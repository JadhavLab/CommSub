 function [SimAnal, ...
           SimAnalTableComponent, ...
           SimAnalTableOverall, ...
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
 

% -----------------------
% Process optional inputs 
% -----------------------
ip = inputParser;
ip.addParameter('bestMatchWithoutReplacement', true);
ip.addParameter('replicateDims', []);  % replicatedims are the dimensions that we replicate our operations over (everything but the last dimensions)
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
    [U,Sig,V] = svd(rr.B_ / rr.V');
    smallestDims = min(size(U,2),size(V,2));
    uNew = nan(cellcount.(source), cellcount.(source));
    vNew = nan(cellcount.(target), cellcount.(target));
    sigNew = nan(cellcount.(source), cellcount.(target));
    %uvNew = zeros(size(uNew,1)+size(vNew,1), smallestDims);
    uNew(sourceInds,1:size(U,2)) = U;
    vNew(targetInds,1:size(V,2)) = V;
    sigNew = Sig;
    % Smallest  number of dims
    Patterns(i).rrSVD = struct(...
        'U',uNew,...
        'V',vNew,...
        'Sig',sigNew,...
        'SVD', {{U,Sig,V}},...
        'UV', [uNew(:,1:smallestDims);vNew(:,1:smallestDims)]);
end

% TODO need to ensure cells mapped into same rows/cols across patterns to compare

% -------------------------------------------------------------------
% Now let's take all partitions and compute a metric, the correlation
% -------------------------------------------------------------------
indices = nd.indicesMatrixForm(Patterns);
if isempty(Opt.replicateDims) && size(indices,2) == 3
    Opt.replicateDims = [1,2];
elseif isempty(Opt.replicateDims) && size(indices,2) == 3
    Opt.replicateDims = [2,3];
end
indRepl = num2cell(indices(:, Opt.replicateDims), 1); % These are the dimensions that we are replicating over
groups  = findgroups(indRepl{:});
uGroups = unique(groups);
SimAnal = struct();
SimAnalTableComponent = table();
SimAnalTableOverall = table();

for g = progress(uGroups','Title','groups')

    releventIndices = indices(g==groups,:);
    R = size(releventIndices,1);

    % -----------------------------------------
    % Create a matrix of the relevent UV groups
    % -----------------------------------------
    for field = progress(["U","V","UV"],'Title','fields')
        for ri = 1:R
            for rj = ri+1:R

                % -------------------------
                % What data will we work on
                % -------------------------
                RI = num2cell(releventIndices(ri,:));
                S1 = Patterns(RI{:});
                RJ = num2cell(releventIndices(rj,:));
                S2 = Patterns(RJ{:});
                S1 = S1.rrSVD;
                S2 = S2.rrSVD;
                tableScaffold = struct(); % For later
                SimAnal(g,ri,rj).partition = releventIndices(1,1);
                SimAnal(g,ri,rj).directionalityNum = releventIndices(1,2);
                for copyField = ["directionality"]
                    SimAnal(g,ri,rj).(copyField) = Patterns(RI{:}).(copyField);
                    tableScaffold.(copyField) = Patterns(RI{:}).(copyField);
                end
                SimAnal(g,ri,rj).pattern1 = Patterns(RI{:}).name;
                SimAnal(g,ri,rj).pattern2 = Patterns(RJ{:}).name;
                tableScaffold.pattern1 = Patterns(RI{:}).name;
                tableScaffold.pattern2 = Patterns(RJ{:}).name;
                tableScaffold.field = field;

                % -----------------------
                % Get similarity of field
                % -----------------------
                principalComponents1 = S1.(field);
                principalComponents2 = S2.(field);
                principalComponents1 = fillmissing(principalComponents1,'constant',0);
                principalComponents2 = fillmissing(principalComponents2,'constant',0);
                product = principalComponents1' * principalComponents2;
                product(all(product == 0, 1),:) = nan;
                product(:,all(product == 0, 2)) = nan;
                SimAnal(g,ri,rj).(field+"corr") = product;


                % ---------------------
                % Iterative min, argmin
                % ---------------------
                P = product;
                pRow = [1:size(product,2)]';
                bestCorr     = nan(1,size(product,2));
                bestCorrPair = nan(1,size(product,2));
                col = 1;
                previous = zeros(size(product,1),1);
                while size(P,1) > 0 && col <= size(P,1)
                    [M,argM] = nanmax(P(:,col),[],1);
                    bestCorr(col) = M;
                    entry= pRow(find((P(:,col) == M) & ~previous,1,'first'));
                    if ~isempty(entry)
                        bestCorrPair(col) = argM;
                    end
                    if Opt.bestMatchWithoutReplacement && ~isnan(M)
                        previous(argM) = 1;
                    end
                    col = col+1;
                end
                SimAnal(g,ri,rj).(field+"bestS")     = bestCorr;
                SimAnal(g,ri,rj).(field+"bestSP") = bestCorrPair;

                % ---------------------
                % Iterative min, argmin
                % ---------------------
                P = product;
                pRow = [1:size(product,2)]';
                worstCorr     = nan(1,size(product,2));
                worstCorrPair = nan(1,size(product,2));
                col = 1;
                previous = zeros(size(product,1),1);
                while size(P,1) > 0 && col <= size(P,1)
                    [M,argM] = nanmin(P(:,col),[],1);
                    worstCorr(col) = M;
                    entry= pRow(find((P(:,col) == M) & ~previous,1,'first'));
                    if ~isempty(entry)
                        worstCorrPair(col) = argM;
                    end
                    if Opt.bestMatchWithoutReplacement && ~isnan(M)
                        previous(argM) = 1;
                    end
                    col = col+1;
                end
                SimAnal(g,ri,rj).(field+"bestD")     = worstCorr;
                SimAnal(g,ri,rj).(field+"bestDP") = worstCorrPair;

                % ----------------------------
                % Reorder UV,U,V by similarity
                % ----------------------------
                order =  SimAnal(g,ri,rj).(field+"bestSP");
                order = order(~isnan(order));
                SimAnal(g, ri, rj).("maxsort"+field) = SimAnal(g,ri,rj).(field+"corr")(order,:);
                % -------------------------------
                % Reorder UV,U,V by dissimilarity
                % -------------------------------
                order =  SimAnal(g,ri,rj).(field+"bestDP");
                order = order(~isnan(order));
                SimAnal(g, ri, rj).("minsort"+field) = SimAnal(g,ri,rj).(field+"corr")(order, :);

                % ----------------------------
                % Create SimAnalTableComponent
                % ----------------------------
                t = tableScaffold;
                t.corrMean = mean(product,'all');
                t.corrMedian = median(product,'all');
                %---------------------------------------------
                t.corrMaxMean = nanmean(bestCorr,'all');
                t.corrMaxMedian = nanmedian(bestCorr,'all');
                %---------------------------------------------
                t.corrMinMean = nanmean(worstCorr,'all');
                t.corrMinMedian = nanmedian(worstCorr,'all');
                %---------------------------------------------
                SimAnalTableOverall = [SimAnalTableOverall; struct2table(t)];
                
                % Create SimAnalTableComponent
                % ----------------------------
                for col = 1:numel(bestCorr)
                    t = tableScaffold;
                    t.component = col;
                    t.corrMax = bestCorr(col);
                    t.corrMin = worstCorr(col);
                    SimAnalTableComponent = [SimAnalTableComponent; struct2table(t)];
                end
                
            end
        end
    end
end
