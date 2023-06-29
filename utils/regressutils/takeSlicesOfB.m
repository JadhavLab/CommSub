function [B_slices] = takeSlicesOfB(extendedB, optDim, withIntercept)

% take the extended mapping B matrix and break it into an cell of matrices
% corresponding to the possible ranks upper bounded by the optimal
% dimensions (with or without intercept)
% if with intercept the rank will be one more than the optDim

B_slices = cell(1,optDim);
numTarget = size(extendedB, 2)/optDim;
numSource = size(extendedB,1)-1;

for i = 1:optDim
    if withIntercept
        B_slices{i} = extendedB(:, (i-1)*numTarget+1:i*numTarget);
    else
        B_slices{i} = extendedB(2:numSource+1, (i-1)*numTarget+1:i*numTarget);
    end
end

end

