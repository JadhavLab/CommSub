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
    cca = cca(i).match;
else
    cca = cca.match;
end

match.activities = cca.activities';
match.activities_source = cca.activities_source';
match.activities_target = cca.activities_target';

match.activities        = match.activities ./ max(match.activities);
match.activities_source = match.activities_source ./ max(match.activities_source);
match.activities_target = match.activities_target ./ max(match.activities_target);

% Constrain each to N samples
match.activities = match.activities(1:N, :);
match.activities_source = match.activities_source(1:N, :);
match.activities_target = match.activities_target(1:N, :);

figure;
tiledlayout(3,1)
nexttile
plot(match.activities);
title('Activities');
nexttile
plot(match.activities_source);
title('Source');
nexttile
plot(match.activities_target);
title('Target');
linkaxes(findall(gcf,'type','axes'),'x');

end

