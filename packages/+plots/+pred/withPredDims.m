%% this script will calculate and plot the preformance as predicitive dimensions
% increases

%% calculate/plot
tab = table();
for i = progress(1:numel(Patterns), 'Title', 'making predDim table')
    P = Patterns(i);
    curr_cvLoss = P.rankRegress.cvLoss;
    curr_rrDim  = P.rankRegress.optDimReducedRankRegress;
    direction   = P.directionality;
    iDataset    = P.iDataset;
    iPartition  = P.iPartition;
    genH        = P.genH_name;
    animal      = P.animal;
    name        = P.name;
    numDimsUsedForPrediction = size(curr_cvLoss,2);
    full_model = 1 - curr_cvLoss(1,end);
    mea        = 1 - curr_cvLoss(1,:);
    err        = curr_cvLoss(2,:);
    dims = 1:numDimsUsedForPrediction;
    mea  = mea(:);
    err  = err(:);
    dims = dims(:);
    optDim   = repmat(curr_rrDim,numDimsUsedForPrediction,1);
    iDataset = repmat(iDataset,numDimsUsedForPrediction,1);
    iP         = repmat(iPartition,numDimsUsedForPrediction,1);
    genH       = repmat(genH,numDimsUsedForPrediction,1);
    animal     = repmat(animal,numDimsUsedForPrediction,1);
    direction  = repmat(direction,numDimsUsedForPrediction,1);
    full_model = repmat(full_model,numDimsUsedForPrediction,1);
    name       = repmat(name,numDimsUsedForPrediction,1);
    numDimsUsedForPrediction = repmat(numDimsUsedForPrediction,numDimsUsedForPrediction,1);
    newtab = table(iDataset,iP,genH,animal,direction,numDimsUsedForPrediction,dims,mea,err,optDim,full_model,name);
    tab = [tab;newtab];
end
tab.fracOptDim = double(tab.optDim)./tab.numDimsUsedForPrediction;
tab.fracDim = tab.dims./tab.numDimsUsedForPrediction;
writetable(tab, figuredefine("tables", "predDim.csv"))

% figure(701); clf
    % subplot(1,3,i)
    % full_model = plots.pred.plotPredictiveDimensions(numDimsUsedForPrediction,...
    %     curr_cvLoss, "optDim", curr_rrDim, "mode", "rr", "color","blue", "normalized", true);
    % hold on
    % xlim([0,8])
    % if j == 1, ax1 = gca;
    % else,      ax2 = gca;
    % end
    % title(P.name)
    % xlim([0,10])
    % linkaxes([ax1, ax2],'y')
    % hold on
% legend ("hpc-hpc", "hpc-pfc")

% Plots
% 1. dim vs performance, error bars are teh dev
% 2. dim vs optDim
