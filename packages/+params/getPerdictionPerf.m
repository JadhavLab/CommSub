function [full_model_performance, full_model_perf_distribution] = ...
           getPerdictionPerf(Patterns, nTarget, nSource, nPartition)
numUsedForPrediction = min(nTarget,nSource);
% make the averaged version
curr_cvLoss = cell(10,2,3);
curr_rrDim = cell(10,2,3);
% if ndims(Patterns) == 4
%     Patterns = Patterns(1);
% end
for p = 1: nPartition
    for i = 1:3
        for j = 1:2
            curr_cvLoss{p,j,i} = Patterns(p,j,i).rankRegress.cvLoss;
            curr_rrDim{p,j,i} = Patterns(p,j,i).rankRegress.optDimReducedRankRegress;
        end
    end
end

full_model_performance = [];
full_model_perf_distribution = [];


for i = 1:3
    for j = 1:2
        full_model = plotPredictiveDimensions(numUsedForPrediction,...
            curr_cvLoss(:,j,i), "optDim", curr_rrDim(:,j,i),"mode", "rr","do_plot", false);
        for p = 1:nPartition
            full_model_single = plotPredictiveDimensions(numUsedForPrediction,...
                curr_cvLoss(p,j,i), "optDim", curr_rrDim(p,j,i),"mode", "rr", "averaged", false,"do_plot", false);
            
            full_model_perf_distribution = [full_model_perf_distribution, full_model_single];
        end
        full_model_performance = [full_model_performance, full_model];
        
    end
end
end


