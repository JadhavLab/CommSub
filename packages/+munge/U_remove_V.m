function varargout = U_remove_V(U, v, n)
% U_remove_V: remove the projection of U on V
%   U: m x n matrix
%   V: m x n matrix
%   n: the number of columns of V to be used
%   varargout{1}: the projection of U on V
%   varargout{2}: the projection of V on U

    if nargin < 3
        n = size(V, 2);
    end
    Projection  = v(:,1:n) * pinv(v(:,1:n)' * v(:,1:n)) * v(:,1:n)';
    u_remove_v = U - U * Projection;
    varargout{1} = u_remove_v;
    if length(varargout) == 2
        v_remove_u = U_remove_V(v, U);
        varargout{2} = v_remove_u;
    end
