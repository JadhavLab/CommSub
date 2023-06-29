function output = subsample(percent, starting_point, stopping_point, input)
% Subsample contiguous range of times
varargout = {};
if iscolumn(input) || isrow(input)
    output = input(starting_point : stopping_point );
elseif ndims(input) == 2
    output = input(starting_point : stopping_point, :);
elseif ndims(input) == 3
    output = input(starting_point : stopping_point, :, :);
end
