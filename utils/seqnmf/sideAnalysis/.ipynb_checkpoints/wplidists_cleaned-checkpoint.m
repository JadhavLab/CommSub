%function [ax, SUMMARY, signals, filters] = crossCorrelate(data, field)
% CROSSCORRELATE applies cross correlation across the different frequency bands
%
% Optionals
% ---------
sparsity = 10e3;    % How sparse to make the yticklabel
field    = 'wpli';  % Which field to plot the cross-correlation of
%optlistassign(who, varargin{:});

groups = findgroups(data.animalcnt(data.subsample_indices), data.epoch(data.subsample_indices));

%% Get averages of bands
%% ---------------------
% Obtain band averages
filter = sieve.bands(data); % Obtain filtration ranges for all meaningful bands
bandavg = @(x) mean(x,2);   % Function for averaging
for band = sieve.bandset(); 
    signals.(band) = bandavg(data.(field)(:,filter.(band))); 
end
SUMMARY.xcorr = xcorr(struct2array(signals));
SUMMARY.zxcorr = zscore(SUMMARY.xcorr,1);
SUMMARY.nzxcorr = bsxfun(@minus, SUMMARY.zxcorr, mean(SUMMARY.zxcorr,2));

%% Compute x and y tick labwls
%% ------------------
[BSx, BSy] = meshgrid(1:numel(sieve.bandset), 1:numel(sieve.bandset));
BS = [BSx(:), BSy(:)];
P = BS; BS=sieve.bandset;
P = BS(P);
P = P(:,1) + "-" + P(:,2);
SUMMARY.xcorr_x = P;M,N] = size(data.(field)); 
SUMMARY.xcorr_y = [-(M-1):1:(M-1)];
SUMMARY.xcorr_yn = 1:2*M-1;
SUMMARY.xcorr_yt = SUMMARY.xcorr_y / median(diff(D.t));
tmp=string(arrayfun(@(x) sprintf('%1.1f',x), round(SUMMARY.xcorr_yt(1:sparsity:end),2),'UniformOutput',false));

for summary = fieldnames(SUMMARY)'
    summary = string(summary);
    if contains(summary, 'xcorr') && ~contains(summary,'_')
        SUMMARY.(summary + "_kurtosis") = kurtosis(SUMMARY.(summary), 1);
        SUMMARY.(summary + "_skew")     = skewness(SUMMARY.(summary),     1);
        SUMMARY.(summary + "_var")      = var(SUMMARY.(summary),      1);
        SUMMARY.(summary + "_mean")     = mean(SUMMARY.(summary),     1);
    end
end

%% Transform into meaningful scale-free cross-correlation
%% ------------------------------------------------------
fig(field + " Zscore Xcorr")
imagesc(SUMMARY.xcorr);
title('xcorr')
cmocean('balance')
colorbar;
% Create X Ticks
xticks(1:numel(P));
xticklabels(SUMMARY.xcorr_x);
xtickangle(45);
% Create Y Ticks
yticks(SUMMARY.xcorr_yn(1:sparsity:end))
yticklabels(SUMMARY.xcorr_yt(1:sparsity:end))
clim([-1, double(quantile(SUMMARY.xcorr(:), 0.99))]);
% Setup horizontal line
h=hline(SUMMARY.xcorr_yn(find(SUMMARY.xcorr_yt==0))); h.Color='black'; h.LineWidth=3;

fig(field + " Zscore Xcorr")
imagesc(SUMMARY.zxcorr);colorbar
title('zxcorr')
cmocean('balance')
clim([-2.5 2.5])
% Create X Ticks
xticks(1:numel(P));
xticklabels(SUMMARY.xcorr_x);
xtickangle(45);
% Create Y Ticks
yticks(SUMMARY.xcorr_yn(1:sparsity:end))
yticklabels(SUMMARY.xcorr_yt(1:sparsity:end))
clim([-2, 2]);
% Setup horizontal line
h=hline(SUMMARY.xcorr_yn(find(SUMMARY.xcorr_yt==0))); h.Color='black'; h.LineWidth=3;
set(gca,'ydir','normal')

tit = field + " Normalized Zscore Xcorr";
fig(tit)
imagesc(smoothdata(SUMMARY.nzxcorr));
title(tit);
cmocean('balance'); clim([-2.5 2.5]); colorbar
% Create X Ticks
xticks(1:numel(P));
xticklabels(SUMMARY.xcorr_x);
xtickangle(45);
% Create Y Ticks
yticks(SUMMARY.xcorr_yn(1:sparsity:end))
yticklabels(SUMMARY.xcorr_yt(1:sparsity:end))
clim([-2, 2]);
% Setup horizontal line
h=hline(SUMMARY.xcorr_yn(find(SUMMARY.xcorr_yt==0))); h.Color='black'; h.LineWidth=3;
set(gca,'ydir','normal')

