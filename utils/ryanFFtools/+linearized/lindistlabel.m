function lindistlabel(varargin)
% Labels an x or y axis with linear distance labels

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('axes', gca);
ip.addParameter('axis', 'x');
ip.addParameter('points', {sprintf('center'),sprintf('choice'),sprintf('side')});
ip.addParameter('labelpositions', [0, 0.4, 1]);
ip.addParameter('textabove', "")
ip.addParameter('allAxVerticals', true);
ip.parse(varargin{:});
opt = ip.Results;
opt.points=string(opt.points);

if contains(opt.axis, 'x')
    ticklabel = @xticklabels;
    tick      = @xticks;
    lim       = @xlim;
    axisShortcut = 'XLim';
end

axes(opt.axes);
pointvals = [];
lims = [min(get(opt.axes, axisShortcut)), max(get(opt.axes, axisShortcut))];
epsilon = 0.025;

tickval = @(alpha) alpha*range(lims) + lims(1);

tickvals = [];
for point = string(opt.points)
    switch point
        case string(opt.points(1))
            tickvals = [tickvals tickval(opt.labelpositions(1))];
        case string(opt.points(2))
            tickvals = [tickvals tickval(opt.labelpositions(2))];
        case string(opt.points(3))
            tickvals = [tickvals tickval(opt.labelpositions(3))];
    end
end

if ~isempty(char(opt.textabove))
    points = string(sprintf('%s\n',opt.textabove)) + opt.points;
else
    points = opt.points;
end

axes(opt.axes);
tick(tickvals);
ticklabel(cellstr(points));
xtickangle(15)

if opt.allAxVerticals
    axs=findobj(gcf,'type','axes');
    for ax = axs'
        axes(ax);
        v=vline(tickvals);
        set(v,'Color','White','LineWidth',4.5,'LineStyle',':');
        alpha(v,0.25);
        %v.Color = [v.Color, 0.4];
    end
end
