function [R, p, heatplt] = findTimeFactorBehaviorCorr(sol, factorName, beh, varargin)

ip = inputParser;
ip.addParameter('colormap','bam');
ip.addParameter('constrainSig', 0.05);
ip.parse(varargin{:})
Opt = ip.Results;

% Create matrix
if isstruct(sol)
    F = sol.factors.(factorName) .* mean(sol.factors.t,1);
else
    F = sol;
end
B = table2array(beh);

% Address
address = [ones(1,size(F,2)), 2*ones(1,size(B,2))];
Q = [F, B]; % times by properties

badfeatures = find(all(isnan(Q),1)); % 
address(badfeatures) = [];
Q(:, badfeatures) = [];
badtimes = find(any(isnan(Q),2)); % 
Q(badtimes, :) = [];
assert(~any(any(isnan(Q),2)));
[R, p] = corrcoef(Q);

% Extract factor-behavior portion
R = R(address == 1, address == 2);
p = p(address == 1, address == 2);

% Heat plot of the R values
behnames = string(beh.Properties.VariableNames);
behnames = replace(behnames, '_', newline);
components = 1:size(F,2);

nexttile;
r = R;
if Opt.constrainSig
    r(p > Opt.constrainSig) = 0;
end
heatplt = heatmap(behnames, components, r);
colorbar;
crameri(Opt.colormap, 'PivotValue', 0)
