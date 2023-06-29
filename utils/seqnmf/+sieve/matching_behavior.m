function data_inds = matching_behavior(t_data, t_behavior)
% SIEVE.MATCHING_BEHAVIOR Determine times in the data. struct that match behavioral times from the behavior struct
%function data_inds = matching_behavior(t_data, t_behavior)

% Inds of the data
data_inds = lookup(double(t_behavior), double(t_data))
