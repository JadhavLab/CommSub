function idphi_series = idphi(x, y, interval)
    % Calculate dx and dy
    dx = [0; diff(x)];
    dy = [0; diff(y)];

    % Calculate phi
    phi = atan2(dy, dx);

    % Preallocate idphi_series
    idphi_series = zeros(1, length(x));

    % Vectorized calculation of IdPhi for each point
    for i = 1:length(x)
        start_index = max(1, i - interval);
        end_index = min(length(x), i + interval);
        idphi_series(i) = sum(abs(diff(phi(start_index:end_index))));
    end
    idphi_series = idphi_series / (2 * interval + 1);
    idphi_series = idphi_series(:);
end
