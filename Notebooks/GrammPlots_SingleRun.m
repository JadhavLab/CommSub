
T = {};
for area1 = ["CA1","PFC"]
    for area2 = ["CA1","PFC"]
    t = table.analyses.pairwiseCorr(Spk, "CA1","PFC");
        T{end+1} = t;
    end
end
T = vertcat(T{:});
T.within = T.source_area == T.target_area;
T.highlow = T.pattern > 3;
T.patternType = categorical(ceil(T.pattern/2));
head(T)
% string(fieldnames(T))';

% fields
%     "source_index"    "target_index"    "correlation_value" "overall_corr_value"    "diff"    "source_area"    "target_area" "pattern"    "diff_norm"    "diff_norm2"    "diff_norm3" "within"    "highlow"    "patternType"    "Properties"    "Row" "Variables"


% Create a gramm object
clf
var = "correlation_value";
if contains(var, "norm")
    q = 0.01;
    inds = T.(var) > quantile(T.(var), q) & T.(var) < quantile(T.(var), 1-q);
else
    inds = true(size(T.(var)));
end
g = gramm('x', T.(var), ...
'color', T.highlow, 'linestyle', T.within, 'lightness', T.within,...
'subset', inds);
% Split the data into columns by patternType
g.facet_grid([], T.patternType);
% Add density plot
g.stat_density("function", "pdf", "kernel", "normal");
% Draw the plot
g.axe_property("YLim", [0 0.1]);
g.draw();
