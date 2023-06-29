function [rawLabels, prettyLabels, labelRepelem, axisCenters] = label(data, field, fieldtype, characerizationtype, varargin)
% LABLER.LABEL

ip = inputParser;
ip.addParameter('marker_type','trajectory_progress'); % trajectory_progress | absolute_location
ip.parse(varargin{:})
opt = ip.Results;

%% Constants
freqFields = {'wpli','S1','S2','phi','phi_cos', 'phi_sin','C','Calt'};

%% Preamble
% type
if nargin < 3
    type = split(field, '_');
    if numel(type) == 1
        fieldtype = type{1};
        characerizationtype = '';
    elseif numel(type) > 1
        fieldtype = type{2};
        characerizationtype = type{1};
    else
        error('Could not parse the field for type infrmation')
    end
end

%% Compute the raw labels
% type specific relabeling
switch characerizationtype
case 'WH'
    labeler_dimension = 1;
    total_dimension = 3;
case 'W'
    labeler_dimension = 1;
    total_dimension = 3;
case ''
    labeler_dimension = 2;
    total_dimension = 2;
otherwise
    error(['Unrecognized labeler type = ' characerizationtype])
end

switch fieldtype

case 'velocity'
    round_precision = 2;
case 'trajdist'
    round_precision = 1;
case 'aglindist'
    round_precision = 1;
case 'reward'
    round_precision = 1;
case 'trajbound'
    round_precision = 1;
case freqFields
    round_precision = 1;
otherwise
    error('Unrecognized labeler type')
end

repelemField = [fieldtype '_repelem'];
if isfield(data, repelemField)
    inherent_repelem = data.(repelemField);
else
    inherent_repelem = 1;
end

% Get raw labels
% --------------
hotencoder_example  = data.(fieldtype);
if ismember(fieldtype, freqFields)
    rawLabels   = data.f(1,:);
elseif isequal(fieldtype, 'reward')
    rawLabels = repelem(1:size(data.orig_reward,2), inherent_repelem);
elseif isequal(fieldtype, 'trajbound')
    rawLabels = repelem(1:size(data.orig_trajbound,2), inherent_repelem);
else
    hotencoder_original = data.(['orig_' fieldtype]);
    [rawLabels] = labeler.internal(data.(field), labeler_dimension, total_dimension, ...
        hotencoder_example, hotencoder_original, ...
        inherent_repelem, round_precision);
end

%% Relabel with pretty labels
% type specific relabeling
switch fieldtype
case 'velocity'
    prettyLabels = string(rawLabels) + " cm/s";
case 'trajdist'
    nLabels = numel(rawLabels)/inherent_repelem;
    nRegions = nLabels/2;
    switch opt.marker_type
    case 'trajectory_progress'
        prettyLabels = ["O_p", "I_p"]' + string(nRegions:-1:1);
        prettyLabels = repelem(prettyLabels, 1, inherent_repelem)';
        prettyLabels = "$" + prettyLabels(:) + "$";
    case 'absolute_location'
        prettyLabels = "a"' + [string(1:nRegions)];
    end
case 'aglindist'
    nLabels = numel(rawLabels);
    nRegions = nLabels/inherent_repelem;
    prettyLabels = repelem(["A"] + string(1:nRegions), 1, inherent_repelem);
case 'reward'
    %nLabels = numel(rawLabels);
    %nRew = 3;
    N = (size(data.orig_reward,2) - 1);
    basepattern = ["Error", "Correct"]' + ["_{t=" + -N:N "}"];
    prettyLabels = repelem(basepattern, 1, inherent_repelem);
case 'trajbound'
    %nLabels = numel(rawLabels);
    %nRew = 3;
    N = size(data.orig_reward,2)/2 - 1;
    basepattern = ["Outbound", "Inbound"]' + ["_{t=" + -N:N "}"];
    prettyLabels = repelem(basepattern, 1, inherent_repelem);
case freqFields
    prettyLabels = string(arrayfun(@(x) sprintf(['%2.' num2str(round_precision) 'f hz'], x), rawLabels, 'UniformOutput', false));
otherwise
    error('Unrecognized labeler type')
end
prettyLabels = prettyLabels(:);
rawLabels    = rawLabels(:);

%% And store repelem output in case user wants to
% cut out the excess fat from these labels
labelRepelem = inherent_repelem;

%% Compute the YCENTERS
% if the labeled axis is to be used with yticks or xticks, one might want to be
% able to specify the centroid of the points that represent different
% numbers 
midpoint   = round(inherent_repelem/2);
axislength = numel(rawLabels);
axisCenters   = midpoint:inherent_repelem:axislength;
axisCenters = axisCenters(:);
rawLabels = rawLabels(axisCenters);
prettyLabels = prettyLabels(axisCenters);
