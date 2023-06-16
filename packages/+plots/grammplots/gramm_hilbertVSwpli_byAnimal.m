%%  from spectra vs. from coherence, try keeping hilbert constant, dim
text_size = 8;
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
corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,1) = gramm( 'x', x,...
    'y', y,...
    'color', categorical(highfilteredsubset.animal));
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

x = highfilteredsubset.percMax_rrDim;
y = lowcoherencesubset.percMax_rrDim;
corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,2) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(highfilteredsubset.animal));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(1,2).facet_grid(categorical(highfilteredsubset.control), categorical(highfilteredsubset.patternAbstract))
g(1,2).geom_point('dodge', 0.5, 'alpha', 0.3);
g(1,2).geom_abline('style','k:');
g(1,2).stat_cornerhist(corner_kws{:});
g(1,2).set_point_options('base_size', point_size);
g(1,2).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(1,2).set_names('x', "High Hilbert" + newline + "Pred Dims", ...
    'y', "Low Coh" + newline + "Pred Dims", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();


x = highfilteredsubset.percMax_rrDim;
y = lowfilteredsubset.percMax_rrDim;
corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,3) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(highfilteredsubset.animal));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(1,3).facet_grid(categorical(highfilteredsubset.control), categorical(highfilteredsubset.patternAbstract))
g(1,3).geom_point('dodge', 0.5, 'alpha', 0.3);
g(1,3).geom_abline('style','k:');
g(1,3).stat_cornerhist(corner_kws{:});
g(1,3).set_point_options('base_size', point_size);
g(1,3).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(1,3).set_names('x', "High Hilbert" + newline + "Pred Dims", ...
    'y', "Low Hilbert" + newline + "Pred Dims", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();

x = lowfilteredsubset.percMax_rrDim;
y = lowcoherencesubset.percMax_rrDim;
corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(2,1) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(lowfilteredsubset.animal));
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
corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(2,2) = gramm('x', x,...
    'y', y,...
    'color', categorical(lowfilteredsubset.animal));
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

x = highcoherencesubset.percMax_rrDim;
y = lowcoherencesubset.percMax_rrDim;
corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(2,3) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(lowfilteredsubset.animal));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(2,3).facet_grid(categorical(lowfilteredsubset.control), categorical(lowfilteredsubset.patternAbstract))
g(2,3).geom_point('dodge', 0.5, 'alpha', 0.3);
g(2,3).geom_abline('style','k:');
g(2,3).stat_cornerhist(corner_kws{:});
g(2,3).set_point_options('base_size', point_size);
g(2,3).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(2,3).set_names('x', "High Coh" + newline + "Pred Dims", ...
    'y', "Low Coh" + newline + "Pred Dims", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')

g.draw();
% axis([g.facet_axes_handles],'square')
%%  from spectra vs. from coherence, try keeping hilbert constant, performance
clf
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
corner_kws = {'edges',-0.15:0.01:0.15, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,1) = gramm( 'x', x,...
    'y', y,...
    'color', categorical(highfilteredsubset.animal));
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

x = highfilteredsubset.full_model_performance;
y = lowcoherencesubset.full_model_performance;
corner_kws = {'edges',-0.15:0.01:0.15, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,2) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(highfilteredsubset.animal));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(1,2).facet_grid(categorical(highfilteredsubset.control), categorical(highfilteredsubset.patternAbstract))
g(1,2).geom_point('dodge', 0.5, 'alpha', 0.3);
g(1,2).geom_abline('style','k:');
g(1,2).stat_cornerhist(corner_kws{:});
g(1,2).set_point_options('base_size', point_size);
g(1,2).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(1,2).set_names('x', "High Hilbert" + newline + "Model Perf", ...
    'y', "Low Coh" + newline + "Model Perf", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')

g.draw();


x = highfilteredsubset.full_model_performance;
y = lowfilteredsubset.full_model_performance;
corner_kws = {'edges',-0.15:0.01:0.15, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,3) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(highfilteredsubset.animal));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(1,3).facet_grid(categorical(highfilteredsubset.control), categorical(highfilteredsubset.patternAbstract))
g(1,3).geom_point('dodge', 0.5, 'alpha', 0.3);
g(1,3).geom_abline('style','k:');
g(1,3).stat_cornerhist(corner_kws{:});
g(1,3).set_point_options('base_size', point_size);
g(1,3).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(1,3).set_names('x', "High Hilbert" + newline + "Model Perf", ...
    'y', "Low Hilbert" + newline + "Model Perf", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')

g.draw();


x = lowfilteredsubset.full_model_performance;
y = lowcoherencesubset.full_model_performance;
corner_kws = {'edges',-0.15:0.01:0.15, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(2,1) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(lowfilteredsubset.animal));
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
corner_kws = {'edges',-0.15:0.01:0.15, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(2,2) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(lowfilteredsubset.animal));
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

x = highcoherencesubset.full_model_performance;
y = lowcoherencesubset.full_model_performance;
corner_kws = {'edges',-0.15:0.01:0.15, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(2,3) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(lowfilteredsubset.animal));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(2,3).facet_grid(categorical(lowfilteredsubset.control), categorical(lowfilteredsubset.patternAbstract))
g(2,3).geom_point('dodge', 0.5, 'alpha', 0.3);
g(2,3).geom_abline('style','k:');
g(2,3).stat_cornerhist(corner_kws{:});
g(2,3).set_point_options('base_size', point_size);
g(2,3).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(2,3).set_names('x', "High Coh" + newline + "Model Perf", ...
    'y', "Low Coh" + newline + "Model Perf", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')

g.draw();
% axis([g.facet_axes_handles],'square')