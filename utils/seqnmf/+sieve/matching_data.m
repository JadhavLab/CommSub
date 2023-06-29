function behavior_inds = matching_data(t_behavior, t_data)
% SIEVE.MATCHING_DATA function data_inds = matching_data(t_behavior, t_data)
% Determine times in the data. struct that match behavioral times from the behavior struct

% Inds of the data
data_inds = lookup(double(t_data), double(t_behavior))

