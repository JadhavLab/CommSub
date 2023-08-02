function [U_avg, V_avg] = avgCCAPartitions(X, N, n)
    % X is your input matrix (time x neurons)
    % N is the number of random partitions

    if nargin < 3
        n = min(5, size(X, 2));
    end

    % Initialize U_avg and V_avg
    U_avg = zeros(size(X, 1), n);
    V_avg = zeros(size(X, 1), n);

    for i = 1:N
        % Create a random partition of X
        idx = randperm(size(X, 2));
        X1 = X(:, idx(1:end/2));
        X2 = X(:, idx(end/2+1:end));

        % Apply CCA to X1 and X2
        [A, B, r, U, V] = canoncorr(X1, X2);

        % Update U_avg and V_avg
        % Here we're using the first pair of canonical variates
        U_avg = U_avg + abs(U(:, n));
        V_avg = V_avg + abs(V(:, n));
    end

    % Average U_avg and V_avg
    U_avg = U_avg / N;
    V_avg = V_avg / N;
end

