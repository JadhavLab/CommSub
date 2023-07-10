function plot_behavior(behavior, behavior_list, varargin)
%PLOT_BEHAVIOR Plot behavior data

% Create an input parser object
ip = inputParser();

% Add optional parameters
ip.addOptional('N', 1000, @(x) isnumeric(x) && isscalar(x));
ip.addOptional('fig', [], @(x) isa(x, 'matlab.ui.Figure'));
ip.addOptional('ax', [],  @(x) isa(x, 'matlab.graphics.axis.Axes') || iscell(x));

% Parse the inputs
ip.parse(varargin{:});

% Get the values of the inputs
N = ip.Results.N;
fig = ip.Results.fig;
ax = ip.Results.ax;

if isempty(fig)
    % Initialize figure
    fig = figure;
end

% Handle axes
if ~isempty(ax) && ~iscell(ax)
    ax = {ax};
end

% Create tiledlayout in given figure if no axes provided
if isempty(ax)
    t = tiledlayout(fig, length(behavior_list),1);
end

use_time = ismember('time', behavior.Properties.VariableNames) && ...
     ismember('time', behavior_list);

% Plot each behavior
for b = 1:length(behavior_list)
    behavior_name = behavior_list{b};
    
    % Check if behavior_name is a valid column in behavior
    if ismember(behavior_name, behavior.Properties.VariableNames)
        % Extract and normalize behavior data
        behavior_data = behavior.(behavior_name);
        behavior_data = behavior_data ./ max(behavior_data);

        % Constrain to N samples
        behavior_data = behavior_data(1:N);
        if use_time
            time = behavior.time(1:N);
        else
            time = 1:N;
        end

        % Select subplot
        if ~isempty(ax) && b <= length(ax)
            subplot(ax{b});
        else
            nexttile(t);
        end

        hold on;
        plot(time, behavior_data)
        title(behavior_name);
    else
        error(['Behavior "' behavior_name '" not found in behavior data.'])
    end
end

% Link axes
linkaxes(findall(gcf,'type','axes'),'x');

end

