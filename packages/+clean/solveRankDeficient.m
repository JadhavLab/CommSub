function A_reconstructed = solveRankDeficient(A, tol)
%  A_reconstructed = solveRankDeficient(A)
%
%  This function takes in a matrix A and returns a matrix A_reconstructed
%  that is the same size as A but with any rank deficient components
%  removed.  This is done by performing a singular value decomposition and
%  removing any singular values that are close to zero.
%
%  Inputs:
%    A  -  A matrix of size M x N

if nargin < 2
    tol = 1e-10;
end

% Perform singular value decomposition
[U, S, V] = svd(A);

% Find any singular values close to zero (you can adjust the tolerance as needed)
zero_sv = abs(diag(S)) < tol;

% Remove the corresponding columns from U, rows from S, and columns from V
U(:, zero_sv) = [];
S(zero_sv, :) = [];
S(:, zero_sv) = [];
V(:, zero_sv) = [];

% Reconstruct the matrix without the rank deficient components
A_reconstructed = U * S * V';

