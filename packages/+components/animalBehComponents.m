function [behavior_running_subspaces] = animalBehComponents(unique_times,...
        throwout_times, total_subspaces, sessionTypePerBin)

running_subspaces = total_subspaces(:,sessionTypePerBin == 1);
%% get the basic behavior info and the subspace components during running times
behavior_running_subspaces = running_subspaces(:,~throwout_times);
behavior_running_subspaces = behavior_running_subspaces (:,unique_times);
