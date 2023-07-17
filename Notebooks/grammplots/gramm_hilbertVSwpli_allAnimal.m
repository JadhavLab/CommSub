%%  from spectra vs. from coherence, try keeping hilbert constant, dim
perc_dim_skewness = zeros(2,3);
text_size = 8;
point_size = 1.5;

clf
fig('hilbert - coherence comparison, perc dimension')
clear g

highfilteredsubset = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-pfc" & T.control == "pattern activity";
lowfilteredsubset = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-pfc" & T.control == "control";
highcoherencesubset = T.generateH == 'fromWpli  fromRipTimes'  ...
    & T.directionality == "hpc-pfc" & T.control == "pattern activity";
lowcoherencesubset = T.generateH == 'fromWpli  fromRipTimes'  ...
    & T.directionality == "hpc-pfc" & T.control == "control";


highfilteredsubset = T(highfilteredsubset,:);
lowfilteredsubset = T(lowfilteredsubset,:);

highcoherencesubset = T(highcoherencesubset,:);
lowcoherencesubset = T(lowcoherencesubset,:);


x = highfilteredsubset.percMax_rrDim;
y = highcoherencesubset.percMax_rrDim;
perc_dim_skewness(1,1) = skewness(y-x);
corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,1) = gramm( 'x', x,...
    'y', y,...
    'color', categorical(highfilteredsubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(1,1).facet_grid(categorical(highfilteredsubset.control), categorical(highfilteredsubset.patternAbstract))
g(1,1).geom_point('dodge', 0.5, 'alpha', 0.3);
g(1,1).geom_abline('style','k:');
g(1,1).stat_cornerhist(corner_kws{:});
g(1,1).set_point_options('base_size', point_size);
g(1,1).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(1,1).set_names('x', "High Hilbert" + newline + "Pred Dims", ...
    'y', "High Coh" + newline + "Pred Dims", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();

y = highfilteredsubset.percMax_rrDim;
x = lowcoherencesubset.percMax_rrDim;
perc_dim_skewness(1,2) = skewness(y-x);
corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,2) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(highfilteredsubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(1,2).facet_grid(categorical(highfilteredsubset.control), categorical(highfilteredsubset.patternAbstract))
g(1,2).geom_point('dodge', 0.5, 'alpha', 0.3);
g(1,2).geom_abline('style','k:');
g(1,2).stat_cornerhist(corner_kws{:});
g(1,2).set_point_options('base_size', point_size);
g(1,2).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(1,2).set_names('x', "Low Coh" + newline + "Pred Dims", ...
    'y', "High Hilbert" + newline + "Pred Dims", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();


y = highfilteredsubset.percMax_rrDim;
x = lowfilteredsubset.percMax_rrDim;
perc_dim_skewness(1,3) = skewness(y-x);
corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,3) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(highfilteredsubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(1,3).facet_grid(categorical(highfilteredsubset.control), categorical(highfilteredsubset.patternAbstract))
g(1,3).geom_point('dodge', 0.5, 'alpha', 0.3);
g(1,3).geom_abline('style','k:');
g(1,3).stat_cornerhist(corner_kws{:});
g(1,3).set_point_options('base_size', point_size);
g(1,3).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(1,3).set_names('x', "Low Hilbert" + newline + "Pred Dims", ...
    'y', "High Hilbert" + newline + "Pred Dims", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();

x = lowfilteredsubset.percMax_rrDim;
y = lowcoherencesubset.percMax_rrDim;
perc_dim_skewness(2,1) = skewness(y-x);
corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(2,1) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(lowfilteredsubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(2,1).facet_grid(categorical(lowfilteredsubset.control), categorical(lowfilteredsubset.patternAbstract))
g(2,1).geom_point('dodge', 0.5, 'alpha', 0.3);
g(2,1).geom_abline('style','k:');
g(2,1).stat_cornerhist(corner_kws{:});
g(2,1).set_point_options('base_size', point_size);
g(2,1).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(2,1).set_names('x', "Low Hilbert" + newline + "Pred Dims", ...
    'y', "Low Coh" + newline + "Pred Dims", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();

x = lowfilteredsubset.percMax_rrDim;
y = highcoherencesubset.percMax_rrDim;
perc_dim_skewness(2,2) = skewness(y-x);
corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(2,2) = gramm('x', x,...
    'y', y,...
    'color', categorical(lowfilteredsubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(2,2).facet_grid(categorical(lowfilteredsubset.control), categorical(lowfilteredsubset.patternAbstract))
g(2,2).geom_point('dodge', 0.5, 'alpha', 0.3);
g(2,2).geom_abline('style','k:');
g(2,2).stat_cornerhist(corner_kws{:});
g(2,2).set_point_options('base_size', point_size);
g(2,2).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(2,2).set_names('x', "Low Hilbert" + newline + "Pred Dims", ...
    'y', "High Coh" + newline + "Pred Dims", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();

y = highcoherencesubset.percMax_rrDim;
x = lowcoherencesubset.percMax_rrDim;
perc_dim_skewness(2,3) = skewness(y-x);
corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(2,3) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(lowfilteredsubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(2,3).facet_grid(categorical(lowfilteredsubset.control), categorical(lowfilteredsubset.patternAbstract))
g(2,3).geom_point('dodge', 0.5, 'alpha', 0.3);
g(2,3).geom_abline('style','k:');
g(2,3).stat_cornerhist(corner_kws{:});
g(2,3).set_point_options('base_size', point_size);
g(2,3).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(2,3).set_names('x', "Low Coh" + newline + "Pred Dims", ...
    'y', "High Coh" + newline + "Pred Dims", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();
axis([g.facet_axes_handles],'square')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')

%%
figure
y = highcoherencesubset.percMax_rrDim;
x = lowcoherencesubset.percMax_rrDim;
perc_dim_skewness(2,3) = skewness(y-x);
corner_kws = {'edges',-0.5:0.025:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(lowfilteredsubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g.facet_grid(categorical(lowfilteredsubset.control), categorical(lowfilteredsubset.patternAbstract))
g.geom_point('dodge', 0.5, 'alpha', 0.3);
g.geom_abline('style','k:');
g.stat_cornerhist(corner_kws{:});
g.set_point_options('base_size', 2);
g.set_text_options('label_scaling', 1.2, 'base_size', text_size);

g.set_names('x', "Low Coh" + newline + "Pred Dims", ...
    'y', "High Coh" + newline + "Pred Dims", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
% g.axe_property('XTickLabelRotation',35)   
% g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();
axis([g.facet_axes_handles],'square')


%%
fig('perc dim difference')
imagesc(perc_dim_skewness)
crameri vik
colorbar
ax = gca;
ax.CLim = [-1 1];
text(0.5,1,'high coherence - high filter ')
text(1.5,1,'low coherence - high filter')
text(2.5,1,'high filter - low filter')

text(0.5,2,'low coherence - low filter')
text(1.5,2,'high coherence - low filter')
text(2.5,2,'high coherence - low coherence')

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

%%  from spectra vs. from coherence, try keeping hilbert constant, performance
clf
perf_skewness = zeros(2,3);

fig('hilbert - coherence comparison, model performance')
clear g

highfilteredsubset = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-pfc" & T.control == "pattern activity";
lowfilteredsubset = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-pfc" & T.control == "control";
highcoherencesubset = T.generateH == 'fromWpli  fromRipTimes'  ...
    & T.directionality == "hpc-pfc" & T.control == "pattern activity";
lowcoherencesubset = T.generateH == 'fromWpli  fromRipTimes'  ...
    & T.directionality == "hpc-pfc" & T.control == "control";


highfilteredsubset = T(highfilteredsubset,:);
lowfilteredsubset = T(lowfilteredsubset,:);

highcoherencesubset = T(highcoherencesubset,:);
lowcoherencesubset = T(lowcoherencesubset,:);


x = highfilteredsubset.full_model_performance;
y = highcoherencesubset.full_model_performance;
perf_skewness(1,1) = skewness(y-x);

corner_kws = {'edges',-0.15:0.01:0.15, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,1) = gramm( 'x', x,...
    'y', y,...
    'color', categorical(highfilteredsubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(1,1).facet_grid(categorical(highfilteredsubset.control), categorical(highfilteredsubset.patternAbstract))
g(1,1).geom_point('dodge', 0.5, 'alpha', 0.3);
g(1,1).geom_abline('style','k:');
g(1,1).stat_cornerhist(corner_kws{:});
g(1,1).set_point_options('base_size', point_size);
g(1,1).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(1,1).set_names('x', "High Hilbert" + newline + "Model Perf", ...
    'y', "High Coh" + newline + "Model Perf", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();

y = highfilteredsubset.full_model_performance;
x = lowcoherencesubset.full_model_performance;
perf_skewness(1,2) = skewness(y-x);

corner_kws = {'edges',-0.15:0.01:0.15, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,2) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(highfilteredsubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(1,2).facet_grid(categorical(highfilteredsubset.control), categorical(highfilteredsubset.patternAbstract))
g(1,2).geom_point('dodge', 0.5, 'alpha', 0.3);
g(1,2).geom_abline('style','k:');
g(1,2).stat_cornerhist(corner_kws{:});
g(1,2).set_point_options('base_size', point_size);
g(1,2).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(1,2).set_names('x', "Low Coh" + newline + "Model Perf", ...
    'y', "High Hilbert" + newline + "Model Perf", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();


y = highfilteredsubset.full_model_performance;
x = lowfilteredsubset.full_model_performance;
perf_skewness(1,3) = skewness(y-x);

corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,3) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(highfilteredsubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(1,3).facet_grid(categorical(highfilteredsubset.control), categorical(highfilteredsubset.patternAbstract))
g(1,3).geom_point('dodge', 0.5, 'alpha', 0.3);
g(1,3).geom_abline('style','k:');
g(1,3).stat_cornerhist(corner_kws{:});
g(1,3).set_point_options('base_size', point_size);
g(1,3).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(1,3).set_names('x', "Low Hilbert" + newline + "Model Perf", ...
    'y', "High Hilbert" + newline + "Model Perf", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();


x = lowfilteredsubset.full_model_performance;
y = lowcoherencesubset.full_model_performance;
perf_skewness(2,1) = skewness(y-x);

corner_kws = {'edges',-0.15:0.01:0.15, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(2,1) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(lowfilteredsubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(2,1).facet_grid(categorical(lowfilteredsubset.control), categorical(lowfilteredsubset.patternAbstract))
g(2,1).geom_point('dodge', 0.5, 'alpha', 0.3);
g(2,1).geom_abline('style','k:');
g(2,1).stat_cornerhist(corner_kws{:});
g(2,1).set_point_options('base_size', point_size);
g(2,1).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(2,1).set_names('x', "Low Hilbert" + newline + "Model Perf", ...
    'y', "Low Coh" + newline + "Model Perf", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();

x = lowfilteredsubset.full_model_performance;
y = highcoherencesubset.full_model_performance;
perf_skewness(2,2) = skewness(y-x);

corner_kws = {'edges',-0.15:0.01:0.15, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(2,2) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(lowfilteredsubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(2,2).facet_grid(categorical(lowfilteredsubset.control), categorical(lowfilteredsubset.patternAbstract))
g(2,2).geom_point('dodge', 0.5, 'alpha', 0.3);
g(2,2).geom_abline('style','k:');
g(2,2).stat_cornerhist(corner_kws{:});
g(2,2).set_point_options('base_size', point_size);
g(2,2).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(2,2).set_names('x', "Low Hilbert" + newline + "Model Perf", ...
    'y', "High Coh" + newline + "Model Perf", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')

g.draw();

y = highcoherencesubset.full_model_performance;
x = lowcoherencesubset.full_model_performance;
perf_skewness(2,3) = skewness(y-x);

corner_kws = {'edges',-0.15:0.01:0.15, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(2,3) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(lowfilteredsubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(2,3).facet_grid(categorical(lowfilteredsubset.control), categorical(lowfilteredsubset.patternAbstract))
g(2,3).geom_point('dodge', 0.5, 'alpha', 0.3);
g(2,3).geom_abline('style','k:');
g(2,3).stat_cornerhist(corner_kws{:});
g(2,3).set_point_options('base_size', point_size);
g(2,3).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(2,3).set_names('x', "Low Coh" + newline + "Model Perf", ...
    'y', "High Coh" + newline + "Model Perf", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();
axis([g.facet_axes_handles],'square')

%%
fig('perf difference')
imagesc(perf_skewness)
crameri vik
colorbar
ax = gca;
ax.CLim = [-1 1];
text(0.5,1,'high coherence - high filter ')
text(1.5,1,'low coherence - high filter')
text(2.5,1,'high filter - low filter')

text(0.5,2,'low coherence - low filter')
text(1.5,2,'high coherence - low filter')
text(2.5,2,'high coherence - low coherence')

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
