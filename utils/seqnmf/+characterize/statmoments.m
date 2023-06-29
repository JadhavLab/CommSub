function statmoments(data, varargin)
% Get statistical moments of a bunch of fields
dim = 1 ;
ip = inputParser;
ip.addParameter('Slog', true);
ip.addParameter('fields',["S1","S2","wpli"])
ip.addParameter('extraFields', string())
ip.parse(varargin{:});
opt = ip.Results;
if ~isempty(opt.extraFields)
    opt.fields = [opt.fields, opt.extraFields];
end

nestable(gca);
count = [0 0];
metrics = {@mean, @var, @skewness, @kurtosis};
for field = string(opt.fields)
    if field == ""
        continue
    end
    count(1)=count(1)+1;
    count(2)=0;
    for func = metrics;
        func=func{1};
        count(2)=count(2)+1;
        nestplot(numel(opt.fields),numel(metrics),{count(1),count(2)});
        if opt.Slog && field.startsWith('S')
            D     = log10(func((data.(field)),dim));
            fname = "log(" + field + ")";
        else
            D = func((data.(field)),dim);
            fname =field;
        end
        plot(data.f, D,'k');
        title(fname + " " + char(func));
    end
end

