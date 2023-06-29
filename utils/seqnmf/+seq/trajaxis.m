function trajaxis(instructions, varargin)
% Takes a set of triplet instructions to label an axis

ip = inputParser;
ip.addParameter('axis','y')
ip.parse(varargin{:});
opt = ip.Results;

% Build the labels
xticklabels = {};
for instuction = instructions'

end
