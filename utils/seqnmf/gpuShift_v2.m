function M = gpuShift_v2(M, k)
% shift 2D array on GPU
    %N = zeros(size(M), 'single', 'gpuArray');
    if k(1) >= 0  && k(2) >= 0
        M(k(1)+1:end, k(2)+1:end) = M(1:end-k(1), 1:end-k(2));
        M(1:k(1), :) = 0;
        M(:, 1:k(2)) = 0;
        M(1:k(1), 1:k(2)) = 0;
    elseif k(1) >= 0  && k(2) < 0
        M(k(1)+1:end, 1:end+k(2)) = M(1:end-k(1), 1-k(2):end);
        M(1:k(1), :) = 0; % Blank out lower register of 1st dim
        M(:, k(2)+end+1:end) = 0; % Blank out upper register of 2nd dim
        M(1:k(1), k(2)+end+1:end) = 0; % Blank out upper/lower square of 1st and 2nd dim
    elseif k(1) < 0  && k(2) < 0
        M(1:end+k(1), 1:end+k(2)) = M(1-k(1):end, 1-k(2):end);
        M(k(1)+end+1:end, :) = 0; % Blank out lower register of 2nd dim
        M(:, k(2)+end+1:end) = 0; % Blank out upper register of 2nd dim
        M(k(1)+end+1:end, k(2)+end+1:end) = 0; % Blank out upper/lower square of 1st and 2nd dim
    elseif k(1) < 0  && k(2) >= 0
        M(1:end+k(1), k(2)+1:end) = M(1-k(1):end, 1:end-k(2));
        M(:, 1:k(2)) = 0; % Blank out lower register of 1st dim
        M(k(1)+end+1:end, :) = 0; % Blank out upper register of 2nd dim
        M(k(1)+end+1:end, 1:k(2)) = 0; % Blank out upper/lower square of 1st and 2nd dim
    else
        error('wrong shift values');
    end
end
