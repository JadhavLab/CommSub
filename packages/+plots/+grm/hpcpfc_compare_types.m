%% from spectra vs. from coherence
clf
fig('hpc-pfc, hilbert vs coherence')
clear g
%T.winSize = cell2mat(T.winSize);
%T = sortrows(T, [16, 21]);
filteredsubset = T.generateH == 'fromSpectra  fromRipTimes'... % NOTE: 2023 ry change to fromSpectra
                    & T.directionality == "hpc-pfc";
coherencesubset = T.generateH == 'fromWpli  fromRipTimes'  ...
                    & T.directionality == "hpc-pfc" ;
filteredsubset = T(filteredsubset,:);
coherencesubset = T(coherencesubset,:);
assert(~isempty(filteredsubset))

x = filteredsubset.percMax_rrDim;
y = coherencesubset.percMax_rrDim;
%subset = filteredsubset.directionality == "hpc-hpc" ;
corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
point_size = 3;

g(1,1) = gramm( ...
            'x', x,...
            'y', y,...
            'color', categorical(filteredsubset.patternAbstract),...
            'lightness', categorical(filteredsubset.control));
% Rassert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g.facet_grid(categorical(filteredsubset.control), categorical(filteredsubset.patternAbstract))
g.geom_point('dodge', 0.5, 'alpha', 0.3);
g.geom_abline('style','k:');
g.stat_cornerhist(corner_kws{:});
g.set_point_options('base_size', point_size);
g.set_text_options('label_scaling', 1.2, 'base_size', 14);

g.set_names('x', "Hilbert" + newline + "Predictive Dims", ...
    'y', "Coh" + newline + "Predictive Dims", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')

x = filteredsubset.full_model_performance;
y = coherencesubset.full_model_performance;
%subset = filteredsubset.directionality == "hpc-pfc" ;
corner_kws = {'edges',-0.16:0.01:0.16, 'aspect',0.6, 'location', [0.2], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,2) = gramm( ... 'subset', subset,...
            'x', x,...
            'y', y,...
            'color',     categorical(filteredsubset.patternAbstract),...
            'lightness', categorical(filteredsubset.control));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(1,2).facet_grid(categorical(filteredsubset.control), categorical(filteredsubset.patternAbstract))
g(1,2).geom_point('dodge', 0.5, 'alpha', 0.3);
g(1,2).geom_abline('style','k:');
g(1,2).stat_cornerhist(corner_kws{:});
g(1,2).set_point_options('base_size', point_size);
g(1,2).set_text_options('label_scaling', 1.2, 'base_size', 14);
g(1,2).set_names('x', "Hilbert" + newline + "Model Performance", ...
    'y', "Coh" + newline + "Model Performance", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.draw();


