function out = single_match_rrr(Patterns, Option, Spk, varargin)
% MATCH Match source and target patterns with rank regressor
% or canonical correlation analysis (CCA) components.
%
%   out = match(Patterns, Option, Spk)
%
%   INPUTS
%       Patterns - struct with fields:
%           X_source     - source patterns
%           X_target     - target patterns
%           X_time       - time vector
%           index_source - index of source patterns
%           index_target - index of target patterns
%           rankRegress  - rank regressor analysis struct
%       Option - struct with fields:
%           method           - 'prod' or 'concat'
%           component_method - 'rank' or 'cca'
%           n_components     - number of components to use
%           smoothing        - samples or param to smooth()
%           verbose          - print out info
%    OUTPUTS
%       out - struct with fields:
%           activities       - matched activities
%           smooth_activities - smoothed activities

ip = inputParser();
ip.addParameter('method', 'prod', @(x) ismember(x, {'prod', 'concat'}));
ip.addParameter('n_components', 3, @isscalar);
ip.addParameter('smoothing', 400); % samples or param to smooth()
ip.addParameter('verbose', false, @islogical);
ip.parse(varargin{:});
Opt = ip.Results;

N                     = Opt.n_components;
out.method            = Opt.method;
out.activities        = [];
out.activities_source = [];
out.activities_target = [];
out.smooth_activities = [];
out.time              = [];
source_index = Patterns.index_source;
target_index = Patterns.index_target;
source       = Spk.spikeCountMatrix(source_index,:);
target       = Spk.spikeCountMatrix(target_index,:);
time         = Spk.timeBinMidPoints;

% Pull out rank regressor
rr = Patterns.rankRegress;
V = rr.V;
d=rr.optDimReducedRankRegress;
B_ = rr.B_ * V(:,1:d)*V(:,1:d)';
if isempty(B_)
    warning('No rank regressor found');
    return % no rank regressor
end
[u,~,v] = svd(B_);
if Opt.verbose
    disp("Size source: " + size(source))
    disp("Size u: " + size(u))
    disp("Size target: " + size(target))
    disp("Size v: " + size(v))
end
activition_source = u(:,1:N)'*source;
activition_target = v(:,1:N)'*target;

% Project onto rank regressor
if strcmp(Opt.method, 'concat')
    activities  = activation_source + activation_target;
elseif strcmp(Opt.method, 'prod')
    % Closer to what occurs in SVD and regression
    % requires product of activations > 0 for match
    activation_source = u(:,1:N)'*source;
    activation_target = v(:,1:N)'*target;
    activities = activation_source .* activation_target;
else
    error('Unknown method');
end

out.activities = activities;
out.activities_source = activation_source;
out.activities_target = activation_target;

% Smooth
smooth_activities = zeros(size(activities));
for i = 1:size(activities,1)
    smooth_activities(i,:) = smooth(activities(i,:), Opt.smoothing);
end
out.smooth_activities = smooth_activities;
out.time = time;
