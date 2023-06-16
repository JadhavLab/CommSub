function [optResult, perf] =...
          optmizeOptions(theta_sp, theta_ripple, Patterns, nTarget, nSource, numPartition)
% this function evaluates how ideal an option struct is 
% by variable result = corr (theta - speed) + corr(theta - ripple)
%                      - avgperf(hpc-pfc)
% avgperf is the averge performance of hpc predicting pfc across the three
% pattern rhythms


% distribution of performance across partitions
if ndims(Patterns) == 4
    Patterns = reshape(Patterns,[numPartition,2,6]);
end
[full_model_performance, dist] = params.getPerdictionPerf(Patterns,...
                                nTarget, nSource, numPartition);
figure(233)
histogram(dist)
perf = sum(full_model_performance(2:2:end))/3;


optResult = theta_sp + theta_ripple - perf;

end

