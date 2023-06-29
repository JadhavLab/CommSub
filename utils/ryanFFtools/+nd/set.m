function X = set(X, ind, x)

ind = num2cell(ind);
X(ind{:}) = x;
