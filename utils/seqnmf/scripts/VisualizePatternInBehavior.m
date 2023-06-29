% Area Wise Pattern
%
% (1) Computes area-wise presence of the pattern progression
% (2) Computes a video of the animal with the pattern wise portions
% overlayed
seqnmf_initialize;
% Get seqnmf
% ---------
load('~/Data/deltacoherence/seqnmf_run_10_lambda=9e-04/master_seqnmf.mat')
% Get best pattern for all times
% ------------------------------
[data.bestpat, data.pattern] = estbestpatterntime(data.data, W, H)


% 1 Region - Pattern analysis
% -------------------
% (Borrow from wenbo: possibly see about grabbing from computer?)


% 2 Video analysis
% -------------------
% (Borrow from wenbo: possibly see about grabbing from computer?)
patternAxis = subplot(8,8,1);
lfpAxis     = subplot(8,8,[6 7 8 10 11 12 14 15 16]);
behaviorAxis = subplot(8,8,[6 7 8 10 11 12 14 15 16]);

