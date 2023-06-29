function isdeficient = isRankDeficient(A)
% test for rank deficiency

sz = size(A);
rk = rank(A);
isdeficient = rk < sz(1);
if isdeficient
    disp('Matrix is rank deficient')
    disp("size(A) = " + sz)
    disp("rank(A) = " + rk)
end
