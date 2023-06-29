 function [tab, ...
           Patterns] = svdSimilarity(Patterns, varargin)
% Carries out svd similarity analysis
%
% And adds relevent svd information to the pattern struct.
%
% Takes in patterns struct and attempts to take in the Patterns struct and
% extract singular value matrices from each of the B matrices. This then takes
% each neuron and figures out how often each neuron is a member of each of the
% singular values.  
%
% ------
% Inputs
% ------
% Patters
%
% -------
% Outputs
% -------
% SimAnal : struct
% SimAnalTableComponent : struct
% SimAnalTableOverall : struct
% Patterns : struct
%
% 
 

% -----------------------
% Process optional inputs 
% -----------------------
ip = inputParser;
ip.addParameter('bestMatchWithoutReplacement', true);
ip.addParameter('replicateDims', []);          % replicatedims are the dimensions that we replicate our operations over (everything but the last dimensions)
ip.addParameter('replicateDimNames', []);      % 
ip.addParameter('patternDims', []);            % 
ip.addParameter('patternDimNames', "pattern"); % 
ip.parse(varargin{:});
Opt = ip.Results;

% ----------------------------------------------------
% Compute max cell indices : we will need to remap each of the pattern
% components to the actual cell numbers they represent
% ----------------------------------------------------
index = nd.fieldGet(Patterns, 'index_target');
directionality = string(nd.fieldGetCell(Patterns, 'directionality'));
directionality = directionality.split('-');
uArea = unique(directionality(:));
cellcount = struct();
for area = uArea(:)'
    dfind = directionality == "hpc";
    dfind = squeeze(dfind(:,:,:,2));
    target_indices = nd.fieldGet(Patterns(dfind), 'index_target');
    cellcount.(area) = max(target_indices(:));
end


% ----------------------------------------------
% First let's just add relelvent svd information
% ----------------------------------------------
for i = progress(1:numel(Patterns),'Title','Singular value decomposition')

    directionality = string(Patterns(i).directionality).split('-');
    source = directionality(1);
    target = directionality(2);
    sourceInds = Patterns(i).index_source;
    targetInds = Patterns(i).index_target;
    rr = Patterns(i).rankRegress;

    [U,Sig,V] = svd(rr.B_ / rr.V'); % run svd

    % Create a space where rows are the real neuron labeels
    smallestDims = min(size(U,2),size(V,2)); % which matrix has the smallest column size? (fewest neurons)
    uNew = nan(cellcount.(source), cellcount.(source)); 
    vNew = nan(cellcount.(target), cellcount.(target));
    sigNew = nan(cellcount.(source), cellcount.(target));
    %uvNew = zeros(size(uNew,1)+size(vNew,1), smallestDims);
    
    % Apply a sort
    % ------------
    [~, uSort] = sort(sourceInds);
    [~, vSort] = sort(targetInds);
    U = U(uSort,:);
    V = V(vSort,:);
    sourceInds = sourceInds(uSort);
    targetInds = targetInds(vSort);
    
    % Place the neuron-agnostic matrices into the neuron stable matrices (where row 2 is neuron 2 for every parition)
    % ---------------------------------------------------------------------------------------------------------------
    uNew(sourceInds,1:size(U,2)) = U;
    vNew(targetInds,1:size(V,2)) = V;
    sigNew = Sig;

    if Opt.reorderAndCleanup
        % Whichever is smallest -- that's our size now
        uCols = all(isnan(uNew), 1) | all(uNew == 0, 1);
        uRows = all(isnan(uNew), 2) | all(uNew == 0, 2);
        uNew(:,uCols) = [];
        uNew(uRows,:) = [];
        vCols = all(isnan(vNew), 1) | all(vNew == 0, 1);
        vRows = all(isnan(vNew), 2) | all(vNew == 0, 2);
        vNew(:,vCols) = [];
        vNew(vRows,:) = [];

        % Whichever is smallest -- that's our size now
        newSize = min([size(uNew); size(vNew)], [], 1);
        uNew   = uNew(:,   1:newSize(2));
        vNew   = vNew(:,   1:newSize(2));
        sigNew = sigNew(:, 1:newSize(2));
    end

    % Smallest  number of dims
    Patterns(i).rrSVD = struct(...
        'Uneurons', sourceInds,...
        'Vneurons', targetInds,...
        'UVneurons', [sourceInds, targetInds],...
        'U',   uNew,...
        'V',   vNew,...
        'Sig', sigNew,...
        'SVD', {{U,Sig,V}},...
        'UV', [uNew(:,1:smallestDims);...
               vNew(:,1:smallestDims)]);
end

% TODO need to ensure cells mapped into same rows/cols across patterns to compare
% Do some computations that trigger an error if rows shared for 

% --------------------------------------------
% Default arguments for Opt params, when empty
% --------------------------------------------
indices = nd.indicesMatrixForm(Patterns);
if isempty(Opt.replicateDims) && size(indices,2) == 3
    Opt.replicateDims = [1,2];
elseif isempty(Opt.replicateDims) && size(indices,2) == 4
    Opt.replicateDims = [1,2,3];
end
if isempty(Opt.replicateDimNames) && size(indices,2) == 3
    Opt.replicateDimNames = ["iPartition", "directionality"];
elseif isempty(Opt.replicateDimNames) && size(indices,2) == 4
    Opt.replicateDimNames = ["generateH", "iPartition", "directionality"];
end
if isempty(Opt.patternDims) 
    szPatterns = size(Patterns);
    Opt.patternDims = szPatterns(ismember(szPatterns, Opt.replicateDims));
end

% -------------------------------------------------------------------
% Now let's take all partitions and compute a metric, the correlation
% -------------------------------------------------------------------
indRepl = num2cell(indices(:, Opt.replicateDims), 1); % These are the dimensions that we are replicating our pattern-wise analysis over
groups  = findgroups(indRepl{:});
uGroups = unique(groups);

for field = ["U", "V"]
    tab.(field) = table();
end

%SimAnalTableComponent = table();
%SimAnalTableOverall = table();

for g = progress(uGroups','Title','groups')

    releventIndices = indices(g==groups,:);
    R = size(releventIndices,1);

    % -----------------------------------------
    % Create a matrix of the relevent UV groups
    % -----------------------------------------
    for field = ["U","V"]
        for pattern = 1:R

            % -------------------------
            % What data will we work on
            % -------------------------
            RI = num2cell(releventIndices(pattern,:));
            S = Patterns(RI{:});
            S = S.rrSVD;
            % Key pieces of information
            partition         = releventIndices(1,1);
            directionalityNum = releventIndices(1,2);
            if isfield(Patterns(RI{:}), 'generateH')
                generateH = Patterns(RI{:}).generateH;
                if contains(generateH, "Wpli")
                    generateH = "wpli";
                elseif contains(generateH, "EEG")
                    generateH = "hilbert";
                end
            else
                generateH = "";
            end
            directionality = Patterns(RI{:}).directionality;
            rotations  = S.(field);
            rotations  = fillmissing(rotations, 'constant', nan);
            rotations  = S.(field);
            magnitudes = S.Sig;

            % -------------------------------------
            % Build the svd table per group, u/v/uv
            % -------------------------------------
            X = [];
            for nIdx = 1:size(rotations, 1)
                for component = 1:size(rotations,2)
                    neuron = S.(field + "neurons");
                    neuron = neuron(nIdx);
                    rotation = rotations(nIdx, component);
                    magnitude = diag(magnitudes);
                    magnitude = magnitude(component);
                    magAdjRotation = magnitude * rotation;
                    x = table(generateH, pattern, directionality, partition, component, neuron, rotation, magnitude, magAdjRotation);
                    X = [X;x];
                end
            end
            tab.(field) = [tab.(field); X];
        end
    end


end % GROUP LOOP
