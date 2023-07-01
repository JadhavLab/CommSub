  function T = getPatternTable(Patterns, Option, varargin)
% Takes pattern struct from a TheScript run and generates the corresponding
% table for that Pattern struct


ip = inputParser;
ip.addParameter('abbreviateGenH',false);
ip.addParameter('additionalAssignments', {});
ip.parse(varargin{:});
Opt = ip.Results;

% Option structure?
% ------------------
if nargin > 1 && ~isempty(Option)
    % If user provides input in the form as an Option structure, we should
    % process those into the table
    if isfield(Option,'winSize')
        Option = rmfield(Option, 'winSize');
    end
    if isfield(Option,'generateH') % this is already a field of patterns
        Option = rmfield(Option, 'generateH');
    end
    for field = string(fieldnames(Option))'
        if isempty(Option.(field))
            Option.(field) = string("");
        end
        if length(Option.(field)) > 1
            Option = rmfield(Option, field);
        end
    end

    rowAddition = struct2table(Option);
else
    rowAddition = table();
end


% Parse each pattern and add to the overall struct
% ------------------------------------------------
% Label the first dimension of each struct
%  if ~isfield(Patterns, "iPartition")
%      disp("iPartition not found in Pattern struct. Which dimension is it?") 
%      I = input();
%      Patterns = nd.dimLabel(Patterns, I, ["iPartition"]);
% end
if ~isfield(Patterns, "iPartition")
    warning('iPartition does not exist. Guessing the dimension and adding....')
    [~,partDim] = max(size(Patterns));
    Patterns = nd.dimLabel(Patterns, partDim, "iPartition");
end
numUsedForPrediction = size(Patterns(1,1,1,1), 2);
T = table();
for pattern = progress(Patterns(:)','Title', 'Creating pattern table')
    
   iPartition = pattern.iPartition;
    patternType = pattern.name;
    directionality = pattern.directionality;
    if isempty(directionality)
        disp("directionality not found in Pattern struct. ..");
        continue
    end
    if isfield(pattern,'generateH')
        generateH = string(pattern.generateH);
    else
        generateH = "";
    end
    names = split(directionality, '-');
    source = names(1);
    target = names(2);
    nSource = size(pattern.X_source,1);
    nTarget = size(pattern.X_target,1);
    animal = pattern.animal;
    %Output properties
    rrDim = pattern.rankRegress.optDimReducedRankRegress;
    if isfield(pattern,'factorAnalysis') && ~isempty(pattern.factorAnalysis) && isfield(pattern.factorAnalysis,'optDimFactorRegress')
        faDim = pattern.factorAnalysis.optDimFactorRegress;
        qOpt = pattern.factorAnalysis.qOpt;
    else
        faDim = nan;
        qOpt = nan;
    end
    
    if isempty(rrDim)
        rankRegressDim = nan;
    end
    
    if isfield(pattern, 'epoch')
        epoch = pattern.epoch;
    else
        epoch = "all";
    end
   
    % if no singular warnings, then we should probably set it to false
    if ~isfield('singularWarning', pattern) || isempty(pattern.singularWarning)
        pattern.singularWarning = false;
    end
    singularWarning = pattern.singularWarning;
   
    maxDim = min(nSource,nTarget);
    percMax_rrDim = rrDim/maxDim;
    percMax_faDim = qOpt/maxDim;
    [full_model_performance,pred_by_perf] = plotPredictiveDimensions...
         (numUsedForPrediction,pattern.rankRegress.cvLoss, "do_plot", false, "averaged", false);
        
    row = table(animal,generateH,epoch,iPartition, source, target, patternType, nSource, nTarget, ...
        directionality, rrDim, percMax_rrDim, qOpt,...
        percMax_faDim, full_model_performance, pred_by_perf(2),...
        singularWarning);
    row = [row, rowAddition];

    T = [T; row];
end

if Opt.abbreviateGenH
    T = table.abbreviateGenH(T);
end

if ~isempty(Opt.additionalAssignments)
    for i = 1:size(Opt.additionalAssignments, 1)
        val = Opt.additionalAssignments{i, 2};    
        T.(Opt.additionalAssignments{i, 1}) = ...
                    repmat(val, height(T)/size(val,1), 1);
    end
end

% ------------------------------------------------
% ERROR CHECK
% Check for singular results in factor analysis
% ------------------------------------------------
factorAnalysisChecksum = sum(T.singularWarning) ...
                         > 0.7 * size(T,1);
if factorAnalysisChecksum
    warning("too many singular results in factor analysis")
    choice = input ("proceed to store results?");
    if choice == "no"
        error ("unsuccessful run")
    end
end
% ------------------------------------------------
