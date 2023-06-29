function new_matrix = unmap_traj(matrix, repelem)
% Undoes the repelem transform

if ismatrix(matrix)
    new_matrix = zeros(size(matrix,1)/repelem, size(matrix,2));
    for i = 1:repelem:size(matrix,1)
        new_matrix(i, :) = mean(matrix(i:i+repelem-1,:),1);
    end
elseif ndims(matrix) == 3
    new_matrix = zeros(size(matrix,1)/repelem, size(matrix,2), size(matrix,3));
    for i = 1:repelem:size(matrix,1)
        new_matrix((i-1)/repelem+1, :, :) = mean(matrix(i:i+repelem-1,:,:),1);
    end
end
