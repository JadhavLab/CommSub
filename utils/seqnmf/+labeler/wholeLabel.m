function [axP, axR, axisCenters, fieldPts] = wholeLabel(data, fields, varargin)
% WHOLELABEL Takes a set of fields and generates a label for that entire spectrum of data
%
% Output
% ------
% axP : Pretty label per point
% axR : Raw label per point
% axisCenters : The points that correspond to the centers of the labels
% fieldPts : The field names per point (would help to hypothetically cut up these large stitched matrices)

% Parse optional input
ip = inputParser;
ip.addParameter('tag', '');
ip.addParameter('removeRepelem', false); % Place each of the computed axes into a single vector
ip.addParameter('collapse', true); % Place each of the computed axes into a single vector
ip.addParameter('downsample', 0);  % Downsample the axes?
ip.addParameter('downsample_specOnly', 0);  % Downsample the axes?
ip.addParameter('protectR', []);   % Preserve R during the the downsample process
ip.addParameter('protectP', string([])); % Preserve P during the downsample process
ip.addParameter('protectField', {}); % Preserve a field from the downsample process
ip.addParameter('excludeFieldFinal', true); % Preserve a field from the downsample process
ip.addParameter('addFieldname', false); % Preserve a field from the downsample process
ip.parse(varargin{:});
opt = ip.Results;

opt.specFields = {'S1','S2', 'wpli', 'C', 'C_alt'};

if opt.tag
    fields = cellfun(@(x) [tag '_' x], fields, 'UniformOutput', false);
end

nFields = numel(fields);

% Construct the axis bounds
[axR, axP, rElem, axisCenters, fieldPts] = deal(cell(1,nFields));
cnt = 0;
for field = fields
    cnt = cnt + 1;
    field = field{1}  ;
    [axR{cnt}, axP{cnt}, rElem{cnt}, axisCenters{cnt}] = labeler.label(data, field);
    assert(isequal(numel(axR{cnt}),numel(axP{cnt}),numel(axisCenters{cnt})))
    fieldPts{cnt} = repmat(string(field), numel(axisCenters{cnt}), 1);
end

nGroups = numel(axP);
axisCenterStarts = arrayfun(@(x) numel(axP{x})*rElem{x}, 1:nGroups);
axisCenterStarts = cumsum(axisCenterStarts);
axisCenterStarts = [0 axisCenterStarts(1:end-1)];

if opt.removeRepelem
    axR         = arrayfun(@(x) axR{x}(1:rElem{x}:end),         1:nGroups, 'UniformOutput', false);
    axP         = arrayfun(@(x) axP{x}(1:rElem{x}:end),         1:nGroups, 'UniformOutput', false);
    axisCenters = arrayfun(@(x) axisCenters{x}(1:rElem{x}:end), 1:nGroups, 'UniformOutput', false);
    fieldPts    = arrayfun(@(x) fieldPts{x}(1:rElem{x}:end),    1:nGroups, 'UniformOutput', false);
end

if opt.downsample || opt.downsample_specOnly % Sparsify the axis points?
    [axPtmp, axRtmp] = deal(axP, axR);
    axP         = arrayfun(@(x) specialDownsample(axP{x},         axPtmp{x}, axRtmp{x}, fields{x}, opt), 1:nGroups, 'UniformOutput', false);
    axR         = arrayfun(@(x) specialDownsample(axR{x},         axPtmp{x}, axRtmp{x}, fields{x}, opt), 1:nGroups, 'UniformOutput', false);
    fieldPts    = arrayfun(@(x) specialDownsample(fieldPts{x},    axPtmp{x}, axRtmp{x}, fields{x}, opt), 1:nGroups, 'UniformOutput', false);
    axisCenters = arrayfun(@(x) specialDownsample(axisCenters{x}, axPtmp{x}, axRtmp{x}, fields{x}, opt), 1:nGroups, 'UniformOutput', false);
end

if opt.addFieldname
    axP = arrayfun(@(x) fields{x} + " " + axP{x}, 1:nFields);
end

if opt.collapse % Collapse points into single vectors? (Will you plot separately or stitched?)
    axP         = cat(1, axP{:});
    axR         = cat(1, axR{:});
    fieldPts    = cat(1, fieldPts{:});
    axisCenters = arrayfun(@(x) axisCenterStarts(x) + axisCenters{x}, 1:nGroups, 'UniformOutput', false);
    axisCenters = cat(1, axisCenters{:});
end

function y = specialDownsample(x, axP, axR, fieldname, opt)
% Handles downsampling with option to protect certain labels from being removed by downsample
if ~ismember(fieldname, opt.protectField)
    y = [];
    if isstring(x)
        y = string(y);
    end
    axP = string(axP);
    if opt.excludeFieldFinal
        switch fieldname
        case opt.specFields
            nX = numel(x) - 1;
        otherwise
            nX = numel(x);
        end
    else
        nX = numel(x);
    end
    for i = 1:nX
        isProtected = ismember(axP(i), opt.protectP) || ismember(axR(i), opt.protectR);
        if opt.downsample_specOnly
            if ismember(fieldname, opt.specFields)
                isDownsample = mod(i-1, opt.downsample_specOnly) == 0;
            else
                isDownsample = 1;
            end
        else
            isDownsample = mod(i-1, opt.downsample) == 0;
        end
        if isDownsample || isProtected;
            y(end+1) = x(i);
        end
    end
    y=y(:);
else
    y = x;
end
