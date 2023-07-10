function to = fieldInit(from, field, varargin)
% to = fieldInit(field, values, varargin)
% see fieldSet for details
% 
% This works like fieldSet, except we're using
% it to initialized a struct from the field
% or list of fields in `field`. The shape of
% `values` (an nd.struct in this case) will
% be used to determine the shape of the output
% struct.

assert(isstruct(from), 'values must be a struct');
to = nd.emptyLike(from);
field = cellstr(field);

for i = progress(1:numel(field), 'Title', 'Initializing..>')
    to = nd.fieldSet(to, field{i}, from, varargin{:});
end

end
