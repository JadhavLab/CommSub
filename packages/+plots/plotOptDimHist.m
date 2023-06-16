function [median_optDim, all_optdim, curr_pattern_dims, median_for_each_pattern, mean_for_each_pattern] = ...
    plotOptDimHist(Patterns, Option, numResult)

all_optdim = [];
median_for_each_pattern = zeros(2,numResult);
mean_for_each_pattern = zeros(2,numResult);

for i = 1:numResult
   
    for j = 1:2
         curr_pattern_dims = [];
        for p = 1:Option.numPartition
            all_optdim = [all_optdim, Patterns(p,j,i).rankRegress.optDimReducedRankRegress];
            curr_pattern_dims = [curr_pattern_dims, Patterns(p,j,i).rankRegress.optDimReducedRankRegress];
        end
        median_for_each_pattern(j,i) = median(curr_pattern_dims);
        mean_for_each_pattern(j,i) = mean(curr_pattern_dims);
    end
    
end

median_optDim = median(all_optdim);

% figure(2333)
% subplot(3,1,1)
% histogram(all_optdim,12);
% 
% subplot(3,1,2)
% plot(median_for_each_pattern)
% subplot(3,1,3)
% plot(mean_for_each_pattern)


