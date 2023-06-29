function [ax, SUMMARY, signals, filters] = crossCorrelatePerformance(data, behavior, field)
% CROSSCORRPERF cross correlate a field with performance

% Optionals
% ---------
sparsity = 10e3;    % How sparse to make the yticklabel
field    = 'wpli';  % Which field to plot the cross-correlation of
optlistassign(who, varargin{:});

behavior

%% Get averages of bands
%% ---------------------
% Obtain band averages
filter = sieve.bands(data); % Obtain filtration ranges for all meaningful bands
bandavg = @(x) mean(x,2);   % Function for averaging
for band = sieve.bandset(); 
    signals.(band) = bandavg(data.(field)(:,filter.(band))); 
end

%% Get perforamnce

SUMMARY.xcorr   = xcorr(struct2array(signals));
SUMMARY.zxcorr  = zscore(SUMMARY.xcorr,1);
SUMMARY.nzxcorr = bsxfun(@minus, SUMMARY.zxcorr, mean(SUMMARY.zxcorr,2));

%% Compute x and y tick labwls
%% ------------------
[BSx, BSy] = meshgrid(1:numel(sieve.bandset), 1:numel(sieve.bandset));
BS = [BSx(:), BSy(:)];
P = BS; BS=sieve.bandset;
P = BS(P);
P = P(:,1) + "-" + P(:,2);
SUMMARY.xcorr_x = P;
tmp=string(arrayfun(@(x) sprintf('%1.1f',x), round(SUMMARY.xcorr_yt(1:sparsity:end),2),'UniformOutput',false));
[M,N] = size(data.(field)); 
SUMMARY.xcorr_y = [-(M-1):1:(M-1)];
SUMMARY.xcorr_yt = SUMMARY.xcorr_y / median(diff(data.t));
SUMMARY.xcorr_yn = 1:2*M-1;
       
