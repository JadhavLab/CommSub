function X = alignToAreaMatrices(B, xinds, yinds, xmax, ymax)
% alignToAreaMatrices(B, xinds, yinds)
%
% take the matrix B and place the xrows at xinds and yrows at yinds

X = zeros(xmax, ymax);
X(xinds, yinds) = B;
