% This script plots the inbound and outbound performance curves of the
% animals, with dashed lines representing epochs. 
% in order to use the script, the behavioral table should be fully produced


fig('perf by epoch/time')
clf
plot(animal_behavior.time, animal_behavior.tperf_inbound, 'Linewidth',2)
hold on
plot(animal_behavior.time, animal_behavior.tperf_outbound, 'Linewidth',2)

groups = findgroups(animal_behavior.epoch);
epoch_changes = find(diff(groups)>0);

xlabel("time (s)"); ylabel("performance");


for i = 1:numel(epoch_changes)
line([animal_behavior.time(epoch_changes(i)) animal_behavior.time(epoch_changes(i))], [0 1], 'color', 'black', 'Linestyle','--')
hold on
legend off
end

legend("inbound", "outbound");