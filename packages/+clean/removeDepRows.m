function [B, jb] = removeDepRows(A, tol)
% Function to remove linearly dependent rows from matrix A

if nargin < 2
    tol = max(size(A)) * eps(norm(A));
end

% Transpose matrix A
AT = A';

% Calculate the reduced row echelon form of the transposed matrix
[RT, jb] = rref(AT);

% Transpose RT back to get R
R = RT';

% Select only the independent rows from the original matrix A to create a new matrix B
B = A(jb, :);
end
