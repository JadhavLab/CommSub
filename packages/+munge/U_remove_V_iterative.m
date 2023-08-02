function varargout = U_remove_V_iterative(U, V, n)

    if nargin < 3
        n = size(V, 2);
    end

    % Initialize output
    u_remove_v = zeros(size(U));
    pinv_VTV = pinv(V(:,1:n)' * V(:,1:n));

    % Loop over rows of U
    for i = 1:size(U, 1)
        u_remove_v(i,:) = U(i,:) - U(i,:) * (V(i,1:n) * (pinv_VTV * V(i,1:n)'));
    end

    varargout{1} = u_remove_v;

    if nargout == 2
        v_remove_u = munge.U_remove_V_iterative(V, U, n);
        varargout{2} = u_remove_v;
    end

end

