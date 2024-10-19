function [Patterns] = rankRegress(Patterns, Option, varargin)
% RANKREGRESS - Reduced rank regression
%   Patterns = rankRegress(Patterns, Option)
%
% INPUTS:
%   Patterns - [1 x 1 struct] see partitionPatterns.m for fields
%   Option   - [1 x 1 struct] see option.setdefaults.m and
%                                 option.defaults.m for all fields
%
% OUTPUTS:
%   Patterns - [1 x 1 struct] see partitionPatterns.m for fields
%                   adds the following fields:
%                   rankRegress - [1 x 1 struct]
%                       cvl - [1 x 1 struct] cross validation loss
%                       cvLoss - [2 x nDims] cross validation loss
%                       optDimReducedRankRegress - [1 x 1 double]
%                       B - [nTarget x nDims] regression weights
%                       B_ - [nSource x nDims] regression weights
%                       V - [nDims x nDims] covariance matrix
%                       muOptLoss - [1 x 1 double] mean optimal loss
%                       stdOptLoss - [1 x 1 double] std optimal loss
%                       muFullLoss - [1 x 1 double] mean full loss
%                       stdFullLoss - [1 x 1 double] std full loss
%                       singlesource_B - [1 x nSource cell] regression
%                           weights for single source prediction
%                       singlesource_optDim - [1 x nSource cell] optimal
%                           dimension for single source prediction
%                       B_rrr - [nSource x nTarget] reduced rank regression
%                           weights
disp("rank regression")
tic

ip = inputParser;
ip.addParameter('verbose', false, @islogical);
ip.addParameter('ploton', false, @islogical);
ip.parse(varargin{:});
verbose = ip.Results.verbose;
ploton = ip.Results.ploton;

if Option.waysOfPartitions ~= 2
    nTarget = size(Patterns(1).X_target,1);
    nSource = min(size(Patterns(1).X_source,1),...
        size(Patterns(1,2,1).X_source,1));
else
    nTarget = unique(arrayfun(@(x) numel(x.index_target), Patterns));
    nSource = unique(arrayfun(@(x) numel(x.index_source), Patterns));
end

% Number of cross validation folds.
numDimsUsedForPrediction = 1:min(nTarget, nSource);
cvNumFolds = double(Option.rankRegress.cvnum);
cvOptions = statset('crossval');
regressMethod = @ReducedRankRegress;
cvOptions.UseParallel = true;

% TODO: LOOK FOR HASH FOLDER -- FIX LOCATION
% RYAN 9/23 - RidgeInit and Scale false may not be good for zscore

if Option.ridgeInit
    disp("RidgeInit == true and Scale == true")
    basefiletitle = store.savename(Option);
    kws = {'RidgeInit', true, 'Scale', true, ...
        'SaveLambda', true, 'ReadLamda', true, ...
        'WRITE_IF_N', 10, 'READ_IF_N', 10, ...
        'verbose', true};
else
    kws = {'RidgeInit', false, 'Scale', false};
end

% ------------------------------------------------------------------
% NOTES
% ------------------------------------------------------------------
% 1. RegressMethod: Used to apply the regression.
% 2. RegressFitAndPredict: Contains two main components:
%    a) Regression function
%    b) Loss evaluation
% 3. RankRegressRoutine:
%    a) Applies the crossval function to RegressFitAndPredict
%    b) Produces the cvLoss matrix and optimal dimension
%    c) Recomputes B, V, and B_ with the optimal dimension
% ------------------------------------------------------------------

disp("Running " + numel(Patterns) + " partitions")
% Initialize cell array to store data for table
tableData = cell(numel(Patterns), 6);
for n = progress(1:numel(Patterns), 'Title', 'RankRegress')
    
    p = Patterns(n);
    
    % when the partition is three-ways, j==1 means same target/source
    % pair and j==2 means diff target/source pair
    % disp("processing rrr for "+p+" partition and the "+i+" pattern
    % "+j+" direction")
    curr_source = double(p.X_source');
    curr_target = double(p.X_target');
    nSource     = size(curr_source, 2);
    nTarget     = size(curr_target, 2);
    
    nan_rows = any(isnan(curr_source), 2) | ...
        any(isnan(curr_target), 2); % setect
    
    if(mean(curr_source, 'all') < 100*eps())
        warning('centering...');
        curr_source = curr_source - mean(curr_source, 2);
        curr_target = curr_target - mean(curr_target, 2);
    end
    patternname = p.directionality + "-" + p.name;
    filetitle = basefiletitle + "-" + patternname;
    if Option.ridgeInit
        extrakws = {'title', filetitle};
    else
        extrakws = {};
    end
    
    cvFun = @(Ytrain, Xtrain, Ytest, Xtest) ...
        RegressFitAndPredict(regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
        numDimsUsedForPrediction, 'LossMeasure', 'NSE',...
        kws{:}, extrakws{:});
    [p.rankRegress.cvl, ...
        p.rankRegress.cvLoss, ...
        p.rankRegress.optDimReducedRankRegress,...
        p.rankRegress.B,...
        p.rankRegress.B_,...
        p.rankRegress.V] ...
        = rankRegressRoutine(cvFun, cvNumFolds, ...
        cvOptions, ...
        curr_target(~nan_rows, :), ...
        curr_source(~nan_rows,:), ...
        numDimsUsedForPrediction);
    
    % (UNCOMMENT TO TEST ONE ITERATION)
    %
    % cvFun(curr_target(~nan_rows, :), ...
    %     curr_source(~nan_rows,:), ...
    %     curr_target(~nan_rows, :), ...
    %     curr_source(~nan_rows,:));
    
    
    % Get last warning and if warning contains rankDeficient, then
    % let's note that
    warn =lastwarn;
    if (contains(warn, 'rankDeficient') || contains(warn, 'rank deficient'))
        p.rankRegress.rankDeficient = true;
    else
        p.rankRegress.rankDeficient = false;
    end
    
    
    % if verbose
    %     disp(newline + "Opt dim: " + p.rankRegress.optDimReducedRankRegress)
    % end
    if ploton
        fig('ploton');
        tiledlayout(1,3);
        nexttile;
        imagesc(curr_source);
        title('Source');
        nexttile;
        imagesc(curr_target);
        title('Target');
        nexttile;
        plot(p.rankRegress.cvLoss(1,:));
        hold on;
        plot(p.rankRegress.optDimReducedRankRegress, ...
            p.rankRegress.cvLoss(1, p.rankRegress.optDimReducedRankRegress), ...
            'r*');
        % Plot shaded area showing standard deviation of loss
        x = 1:size(p.rankRegress.cvLoss,2);
        keyboard
        y = p.rankRegress.cvLoss(1,:);
        e = p.rankRegress.cvLoss(2,:);
        x2 = [x, fliplr(x)];
        inBetween = [(y-e), fliplr((y+e))];
        fill(x2, inBetween, 'b');
        alpha(.1);
        hold off;
        title('Cross validation loss');
        xlabel('Number of dimensions');
        ylabel('Loss');
        drawnow;
    end
    
    p.rankRegress.muOptLoss   = p.rankRegress.cvLoss(1, p.rankRegress.optDimReducedRankRegress);
    p.rankRegress.stdOptLoss  = p.rankRegress.cvLoss(2, p.rankRegress.optDimReducedRankRegress);
    p.rankRegress.muFullLoss  = p.rankRegress.cvLoss(1, end);
    p.rankRegress.stdFullLoss = p.rankRegress.cvLoss(2, end);
    disp("Direction " + p.directionality + " with optDim " + p.rankRegress.optDimReducedRankRegress)
    
    % Collect data for table
    tableData{n, 1} = Option.animal;
    tableData{n, 2} = p.name;
    tableData{n, 3} = p.directionality;
    tableData{n, 4} = p.iPartition;  % iPartition
    tableData{n, 5} = p.rankRegress.optDimReducedRankRegress;
    tableData{n, 6} = p.rankRegress.muOptLoss;
    
    % Calculate percOptDim (optimal dimension / possible dimension)
    possibleDim = min(nSource, nTarget);
    percOptDim = (p.rankRegress.optDimReducedRankRegress / possibleDim) * 100;
    tableData{n, 7} = percOptDim;
    
    % Add r2 value (performance)
    r2 = 1 - p.rankRegress.muOptLoss;  % Assuming NSE loss, r2 = 1 - NSE
    tableData{n, 8} = r2;
    
    % Single neuron prediction
    if Option.analysis.singleNeuronPrediction
        
        B_singleprediction   = cell(1,nSource);
        dim_singleprediction = cell(1,nSource);
        
        for k = 1:nSource
            curr_singlesource = curr_source(:,k);
            if clean.zeroFiring(curr_singlesource)
                continue;
            end
            try
                [~,~, dim_singleprediction{k}, B_singleprediction{k},~,~] = ...
                    rankRegressRoutine(cvFun,...
                    cvNumFolds, ...
                    cvOptions, ...
                    curr_target, ...
                    curr_singlesource,...
                    numDimsUsedForPrediction);
            catch
                dim_singleprediction{k} = nan;
                B_singleprediction{k} = [];
            end
            summaryTable = table(tableData(:,1), tableData(:,2), tableData(:,3), cell2mat(tableData(:,4)), ...
                cell2mat(tableData(:,5)), cell2mat(tableData(:,6)), cell2mat(tableData(:,7)), cell2mat(tableData(:,8)), ...
                'VariableNames', {'Animal', 'Pattern', 'Directionality', 'iPartition', ...
                'OptDim', 'MuOptLoss', 'PercOptDim', 'R2'});
            p.rankRegress.B_rrr = getReducedB_(p.rankRegress.B,...
                p.rankRegress.V, ...
                nSource, ...
                nTarget, ...
                p.rankRegress.optDimReducedRankRegress);
        end
        
        
        Patterns(n) = p;
    end
    
    
end

disp("rank regression done in "+toc+" seconds")

% Create summary table
summaryTable = table(tableData(:,1), tableData(:,2), tableData(:,3), cell2mat(tableData(:,4)), ...
    cell2mat(tableData(:,5)), cell2mat(tableData(:,6)), ...
    'VariableNames', {'Animal', 'Pattern', 'Directionality', 'iPartition', ...
    'OptDim', 'MuOptLoss'});


% Display the first few rows of the summary table
disp(head(summaryTable));

% Optionally, you can save the table to a file
savename = store.savename(Option);
filename = sprintf('rank_regress_%s.csv', savename);
filepath = fullfile(hashdefine(), filename);
writetable(summaryTable, filepath);

end
