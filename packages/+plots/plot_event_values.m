function plot_event_values(out_struct)
%PLOT_EVENT_VALUES Plot the event values
% see cca.event_analysis for details on the out_struct

% Get the event_u_values and event_v_values
event_u_values = out_struct.event_u_values;
event_v_values = out_struct.event_v_values;

% Get the sizes of the 3D matrices
[numRows, numSamples, numCols] = size(event_u_values);

% Create a new figure
figure;

% Create a 3D tiled layout
Nh = numRows;
Nw = numCols;
gap    = [0.01 0.03];
marg_h = [0.1 0.01];
marg_w = [0.01 0.01];
ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w);

% Loop over the rows and columns
for row = 1:numRows
    for col = 1:numCols
        % Select the subplot for this row and column
        t = ha((row - 1) * numCols + col);
        % Set axes
        axes(t);

        % Get the u and v values for this row and column
        u_values = squeeze(event_u_values(row, :, col));
        v_values = squeeze(event_v_values(row, :, col));

        % Plot the u and v values
        scatter(u_values, v_values);
        axis equal;
        % hline at 0
        xline(0);
        % vline at 0
        yline(0);
        alpha(0.25);

        % X=Y line with a text label "communication subspace" along
        % the line (45 degrees, 1.5 units from the origin)
        line([-3 3], [-3 3], 'Color', 'red', 'LineStyle', '--');
        tt=text(-1, -1.3, 'communication subspace', 'Color', 'red');
        tt.Rotation = 45;


        % Set the title
        title(sprintf('Pattern: %d, Component: %d', row, col));
    end
end
linkaxes(ha, 'xy');
set(findobj(gcf,'type','axes'),'XLim',[-3 3],'YLim',[-3 3]);

set(gcf, 'Position', get(0, 'Screensize'));
end
