% GRAMM plots for a single run of data

if ~munge.detectZscore(Spk.spikeRateMatrix) && ...
    ~any(Spk.spikeRateMatrix(:) < 0)
    disp(" Z-scoring ")
    if ~isfield(Spk, 'muFR')
        Spk.muFR  = mean(Spk.spikeRateMatrix, 2);
        Spk.stdFR = std(Spk.spikeRateMatrix, 0, 2);
    end
    Spk.spikeRateMatrix  = zscore(Spk.spikeRateMatrix,  0, 2);
    Spk.spikeCountMatrix = zscore(Spk.spikeCountMatrix, 0, 2);
    Spk.avgFR = mean(Spk.spikeRateMatrix, 2);
end

T = {};
corrs = struct();
for area1 = ["CA1","PFC"]
    for area2 = ["CA1","PFC"]
    [t, corrs.(area1+area2), overall]  = table.analyses.pairwiseCorr(Spk, "CA1","PFC");
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

% ISSUE: GETTING ALMOST SAME CORRELATION FOR ALL i AND j
% Create a gramm object
clf
var = "diff";
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
g.set_stat_options("alpha", 0.5, "nboot", 1000);
g.stat_density("function", "pdf")
g.stat_summary("type", "sem", "geom", "area")
% Draw the plot
% g.axe_property("YLim", [0 0.2]);
g.draw();

% correlations minus overall correlation
figure;
tiledlayout(4,6);
F = fieldnames(corrs);
for i = 1:4
for j = 1:6
    nexttile
    imagesc(corrs.(F{i}){j} - overall)
    title(F{i})
    axis square
    cmocean('balance')
    clim([-0.1 0.1])
end
end
