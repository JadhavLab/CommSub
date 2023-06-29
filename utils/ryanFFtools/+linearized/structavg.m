function T = structavg(S, dim, varargin)
% linearized.structavg applies proper averaging method for different fields across struct dimensions dim

% How will the average be performed
avgfuns = struct();
for field = ["C","wpli","S1","S2","t","f"]
    avgfuns.(field) = @nanmean;
end
for field = ["phi"]
    avgfuns.(field) = @circ_mean;
end

% Perform averaging
T = nd.structfun(S, dim, 'fun', avgfuns, varargin{:});
