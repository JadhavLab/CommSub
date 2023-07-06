function activities(cca, varargin)

% Create an input parser object
ip = inputParser();

% Add optional parameters
ip.addOptional('i', [], @(x) isnumeric(x) && isscalar(x));
ip.addOptional('N', 1000, @(x) isnumeric(x) && isscalar(x));

% Parse the inputs
ip.parse(varargin{:});

% Get the values of the inputs
i = ip.Results.i;
N = ip.Results.N;

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

% Compute multitaper spectrogram
params.Fs = 1/median(diff(match.time)); % sampling frequency
params.tapers = [3 5]; % taper parameters, [TW K]
params.pad = 0; % padding factor, can be changed depending on the data
movingwin = [1.5, 0.05]; % moving window [window_length step_size]

disp('Computing multitaper spectrogram...');
[S_activities,t,f]        = mtspecgramc(match.activities,movingwin,params);
[S_activities_source,~,~] = mtspecgramc(match.activities_source,movingwin,params);
[S_activities_target,~,~] = mtspecgramc(match.activities_target,movingwin,params);
disp('Done.');

S_activities = S_activities(1:N,:,1);
S_activities_source = S_activities_source(1:N,:,1);
S_activities_target = S_activities_target(1:N,:,1);

figure;
tiledlayout(3,1)

% Multitaper spectrum of activities
nexttile
imagesc(t,f,10*log10(S_activities)');
axis xy;
title('Multitaper Spectrogram Activities');
colorbar;

% Multitaper spectrum of activities_source
nexttile
imagesc(t,f,10*log10(S_activities_source)');
axis xy;
title('Multitaper Spectrogram Source');
colorbar;

% Multitaper spectrum of activities_target
nexttile
imagesc(t,f,10*log10(S_activities_target)');
axis xy;
title('Multitaper Spectrogram Target');
colorbar;

linkaxes(findall(gcf,'type','axes'),'x');

end

