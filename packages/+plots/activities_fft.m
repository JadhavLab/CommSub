function activities(cca, varargin)

% Create an input parser object
ip = inputParser();

% Add optional parameters
ip.addOptional('i', [], @(x) isnumeric(x) && isscalar(x));
% ip.addOptional('N', 1000, @(x) isnumeric(x) && isscalar(x));

% Parse the inputs
ip.parse(varargin{:});

% Get the values of the inputs
i = ip.Results.i;

% Check if i is provided
if ~isempty(i)
    match = cca(i).match;
else
    match = cca.match;
end

match.activities = match.activities';
match.activities_source = match.activities_source';
match.activities_target = match.activities_target';

match.activities        = match.activities ./ max(match.activities);
match.activities_source = match.activities_source ./ max(match.activities_source);
match.activities_target = match.activities_target ./ max(match.activities_target);

% Compute multitaper spectrum
params.Fs = 1/median(diff(match.time)); % sampling frequency
params.tapers = [3 5]; % taper parameters, [TW K]
params.pad = 0; % padding factor, can be changed depending on the data

[S_activities,f1] = mtspectrumc(match.activities,params);
[S_activities_source,f2] = mtspectrumc(match.activities_source,params);
[S_activities_target,f3] = mtspectrumc(match.activities_target,params);

% Smooth each
origsz = size(S_activities);
for i = 1:origsz(2)
    S_activities(:,i) = smooth(S_activities(:,i), 100);
    S_activities_source(:,i) = smooth(S_activities_source(:,i), 100);
    S_activities_target(:,i) = smooth(S_activities_target(:,i), 100);
end
newsz = size(S_activities);

figure;
tiledlayout(3,1)

% Multitaper spectrum of activities
nexttile
plot(f1, 10*log10(S_activities));
title('Multitaper Spectrum Activities');

% Multitaper spectrum of activities_source
nexttile
plot(f2, 10*log10(S_activities_source));
title('Multitaper Spectrum Source');

% Multitaper spectrum of activities_target
nexttile
plot(f3, 10*log10(S_activities_target));
title('Multitaper Spectrum Target');

linkaxes(findall(gcf,'type','axes'),'x');

end

