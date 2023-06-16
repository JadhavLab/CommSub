function [learning_times, stable_times, learning_phase, ...
         stable_phase, learning_rate, stable_rate, throwout_times] ...
          = calculateLearning(times, perf, quantile)

      %DEPRECATED: NOT REALLY USED
% this function calculates the learning of the animal by its perfromance.
% learning phase is defined by less than 80% of the best performance
% stable phase is defined by the other way

% rate of learning is calculated by taking the derivative of the
% conseuctive time points

throwout_times = isnan(perf);
perf = perf(~throwout_times);
times = times(~throwout_times);
max_perf = max(perf);
threshold = quantile*max_perf;

learning_times = perf<threshold;
stable_times = perf>=threshold;

learning_phase = perf(learning_times);
stable_phase = perf(stable_times);

all_rate(1) = 0;
for i = 1:numel(perf)-1
    all_rate(end+1) = abs((perf(i+1)-perf(i))/(times(i+1)-times(i)));
end

learning_rate = all_rate(learning_times);
stable_rate = all_rate(stable_times);
end

