function packedRawData(data, varargin)
% PACKEDRAWDATA plots data from the field packed .data struct field

ip = inputParser;
ip.addParameter('sep', false)
ip.addParameter('ax', [])
ip.parse(varargin{:});
opt = ip.Results;

if isempty(opt.ax)
    ax = gca;
else
    ax = opt.gca;
end

if opt.sep
    ax = nestable(ax);
    nFields = numel(data.fields);
    cFields = 0;
    AX = {};
    for field = data.fields
        cFields = cFields + 1;
        AX{cFields} = nestplot(nFields, 1, cFields, 'yspacing', 0);
        if cFields == nFields
            xtcs = (1:T) * 1/data.samprate;
        else
            xticks([])
        end
        [raw, pretty, ~, ycenters] = labeler.label(data, ['W_', field{1}]);
    end
else
end
