function out = single_match_cca(Patterns, Option, Spk, varargin)
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

if isempty(Patterns.cca)
    warning('No canonical correlation analysis found');
    return
end

% Pull out cca
cca = Patterns.cca;
A = cca.a;
B = cca.b;
if isempty(A) || isempty(B)
    warning('No canonical correlation analysis found');
    return
end
    
activation_source = A(:,1:N)'*source;
activation_target = B(:,1:N)'*target;

% Project onto cca
if strcmp(Opt.method, 'concat')
    activities = activation_target + activation_source;
    activities = activities - mean(activities,2);
elseif strcmp(Opt.method, 'prod')
    activities = activation_source .* activation_target;
else
    error('Unknown method');
end

activities = activities - mean(activities,2);
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

end

