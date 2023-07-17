if ~exist('T', 'var')
    PreambFigs; % Load data
end

%% COMPARISON OF FULL MODEL PERF AND PREDICTIVE DIMS
plots.grm.hpcpfc_compare_types(T, "power", "coherence");
plots.grm.hpcpfc_compare_types(T, "power", "wpli");
plots.grm.hpcpfc_compare_types(T, "coherence", "wpli");

%% Average pattern maxDimPerc single direction, all patterns
figure(5)
clf
subset    = T.directionality == "hpc-pfc";
hpcsubset = T.directionality == 'hpc-hpc';
pfcsubset = T.directionality == 'hpc-pfc';

% x = categorical(T(hpcsubset, "rrDim"));
% y = T(pfcsubset, "rrDim");
x = categorical(T.patternType);
y = T.percMax_rrDim; 
% y = T.full_model_performance;
% g = gramm(  'x', x,...
%     'y', y);
% g.geom_point();
% g.stat_cornerhist('edges',-4:0.2:4,'aspect',0.6);
g = gramm('subset', subset, ...
    'x', x ,...
    'y', y,...
    'Color', categorical(T.patternAbstract), ...
    'lightness', categorical(T.control));

g.geom_jitter('alpha', 0.8, 'width', 0.6)
g.stat_summary('type', 'sem','geom','black_errorbar')
%g.stat_summary('type', 'sem','geom','bar')
g.set_point_options('base_size', 20);
g.set_text_options('label_scaling', 1.5, 'base_size', 14);

g.set_names('x', 'Pattern Type', ...
    'y', '#(Predictive Dims)/#(Max. Pred. Dims)', ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.set_title(" FilteredEEG Example Pattern Dimensionalitiy" + newline() + "HPC-PFC");

g.axe_property('XTickLabelRotation',35)
g.draw()

% ------------------------------------------------------------

%%$% Cornerhist hpc versus pfc predictive
hpcsubset =   T.directionality == "hpc-hpc";
pfcsubset =  T.directionality == "hpc-pfc" ;
hpcsubset = T(hpcsubset,:);
pfcsubset = T(pfcsubset,:);
x = hpcsubset.rrDim;
y = pfcsubset.rrDim;

clf
f=figure(5)
corner_kws = {'edges',-10:0.5:10, 'aspect',1, 'location', [16], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
% LOCATION SETS WHERE THE CORNER HISTS ARE PLACED
set(f,'Position',get(0,'ScreenSize'));
g = gramm( ...
            'x', x,...
            'y', y,...
            'subset', hpcsubset.control == "control")
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g.facet_grid(categorical(hpcsubset.patternAbstractSymbol), [])
g.geom_jitter('alpha', 0.01, 'width',0.35, 'height',0.35);
g.stat_cornerhist(corner_kws{:});
g.set_point_options('base_size',5);
g.set_text_options('label_scaling', 1.5, 'base_size', 10);
g.set_names('x', "HPC" +newline+ "predictive dimensions", ...
    'y', "PFC" +newline+ "predictive dimensions", ...
    'Color', 'Pattern', ...
    'row','',...
    'Lightness', 'Treatment/Control')
g=g.set_color_options('chroma',0);
g.set_text_options('interpreter','latex','base_size',10)
%set(g.results.geom_jitter_handle,'MarkerSize',5)
g.update('subset', hpcsubset.control ~= "control",...
         'color', categorical(hpcsubset.patternAbstract));
g.set_color_options();
g.stat_cornerhist(corner_kws{:});
g.geom_abline('style','k:');
% g.axe_property('XTickLabelRotation',35, 'axis','square')
g.geom_jitter('alpha', 0.30, 'width',0.35, 'height',0.35);
g.draw()
%set(g.results.geom_jitter_handle,'MarkerSize',5);

% ------------------------------------------------------------

x = hpcsubset.percMax_rrDim;
y = pfcsubset.percMax_rrDim;
f=figure(1005)
corner_kws = {'edges',-0.8:0.05:0.8, 'aspect',1, 'location', [1], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
% LOCATION SETS WHERE THE CORNER HISTS ARE PLACED
set(f,'Position');
g = gramm( ...
            'x', x,...
            'y', y,...
            'subset', hpcsubset.control == "control")
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g.facet_grid(categorical(hpcsubset.patternAbstractSymbol), [])
g.geom_jitter('alpha', 0.01, 'width',0.02, 'height',0.02);
g.stat_cornerhist(corner_kws{:});
g.set_point_options('base_size', 5);
g.set_text_options('label_scaling', 1.5, 'base_size', 10);
g.set_names('x', "HPC" +newline+ "predictive dimensions", ...
    'y', "PFC" +newline+ "predictive dimensions", ...
    'Color', 'Pattern', ...
    'row','',...
    'Lightness', 'Treatment/Control')

g=g.set_color_options('chroma',0);
g.set_text_options('interpreter','latex','base_size',10)
%set(g.results.geom_jitter_handle,'MarkerSize',5)
g.draw()
results1 = g.results;
g.update('subset', hpcsubset.control ~= "control",...
         'color', categorical(hpcsubset.patternAbstract));
g.set_color_options();
g.stat_cornerhist(corner_kws{:});
g.geom_abline('style','k:');
% g.axe_property('XTickLabelRotation',35, 'axis','square')
g.geom_jitter('alpha', 0.50, 'width',0.35, 'height',0.35);
g.draw()
%set(g.results.geom_jitter_handle,'MarkerSize',5);

% ------------------------------
% Permuation tests
% ------------------------------
[G, patAb] = findgroups(hpcsubset.patternAbstract);
p = zeros(size(unique(G)));
for gg = unique(G)'
    t_pat = hpcsubset(G == gg & hpcsubset.control == "pattern activity", :).percMax_rrDim - pfcsubset(G==gg & pfcsubset.control=="pattern activity",:).percMax_rrDim;
    t_ctrl = hpcsubset(G == gg & hpcsubset.control == "control", :).percMax_rrDim - pfcsubset(G==gg & pfcsubset.control=="control",:).percMax_rrDim;
    p(gg) = permutationTest(t_pat, t_ctrl, 1000, 'plotresult', 1, 'meanfunc', @nanmedian);
    title(patAb(gg));
end

% ------------------------------
% Add median to each corner_hist
% ------------------------------
[G, patAb, ctrl] = findgroups(hpcsubset.patternAbstract, hpcsubset.control);
[~,~,uPatAb] = unique(patAb);
[~,~,uCtrl]  = unique(ctrl);
med_hpc = splitapply(@median, hpcsubset.percMax_rrDim, G); % get median of each group
[G, patAb, ctrl] = findgroups(pfcsubset.patternAbstract, pfcsubset.control);
[~,~,uPatAb] = unique(patAb);
[~,~,uCtrl] = unique(ctrl);
med = splitapply(@median, hpcsubset.percMax_rrDim - pfcsubset.percMax_rrDim, G); % get median of each group
for gg = unique(G)'
    corner_hist_handle = g.results.stat_cornerhist(uPatAb(gg)).child_axe_handle;
    if ctrl(gg) == "control"
        corner_hist_color  = [0.5, 0.5, 0.5];
    else
        corner_hist_color  = g.results.stat_cornerhist(uPatAb(gg)).bar_handle.FaceColor/2;
    end
    Y = ylim(corner_hist_handle);
    l = line(corner_hist_handle, [med(gg) med(gg)], Y*1.1,  'Color', corner_hist_color, 'LineWidth', 3, 'LineStyle', ':'); 
end

% ------------------------------------------------------------


%% Cornerhist hpc versus pfc FA qOpt
clf
figure(7)
hpcsubset = T.directionality == 'hpc-hpc';
pfcsubset = T.directionality == 'hpc-pfc';
hpcsubset = T(hpcsubset,:);
pfcsubset = T(pfcsubset,:);
x = hpcsubset.qOpt;
y = pfcsubset.qOpt;
subset =  ~hpcsubset.singularWarning;
g = gramm(  'subset', subset,...
            'x', x,...
            'y', y,...
            'color', categorical(hpcsubset.patternAbstract),...
            'lightness', categorical(hpcsubset.control));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g.facet_grid(categorical(hpcsubset.control), categorical(hpcsubset.patternAbstract))
% g.geom_point('dodge', 0.5, 'alpha', 0.3, 'jitter',0.1);
g.geom_abline('style','k:');
g.stat_cornerhist('edges',-4:0.5:4, 'aspect',1.2);
g.set_point_options('base_size', 10);
g.set_text_options('label_scaling', 1.5, 'base_size', 14);
g.set_names('x', "HPC" +newline+ "regional dimensions", ...
    'y', 'PFC regional dimensions', ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('XTickLabelRotation',35)
g.draw()


%% Cornerhist dimension complexity vs. regression -- PFC more complex but required fewer prediction dimensions
clf
figure(8)
hpcsubset = T.directionality == "pfc-hpc";
pfcsubset = T.directionality == "pfc-pfc";
hpcsubset = T(hpcsubset,:);
pfcsubset = T(pfcsubset,:);
% x = hpcsubset.qOpt./hpcsubset.rrDim;
x1 = pfcsubset.qOpt;
y1 = pfcsubset.rrDim;
x2 = hpcsubset.qOpt;
y2 = hpcsubset.rrDim;
subset = pfcsubset.generateH == "fromFilteredEEG  fromRipTimes" & ~pfcsubset.singularWarning;
g1 = gramm(  'subset', subset,...
            'x', x1,...
            'y', y1,...
            'color', categorical(hpcsubset.patternAbstract),...
            'lightness', categorical(hpcsubset.control));


% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g1.facet_grid(categorical(hpcsubset.control), categorical(hpcsubset.patternAbstract))
g1.geom_point('dodge', 0.5, 'alpha', 0.3);
g1.geom_abline('style','k:');
g1.stat_cornerhist('edges',-4:0.5:4, 'aspect',1.2);
g1.set_point_options('base_size', 10);
g1.set_text_options('label_scaling', 1.5, 'base_size', 14);

g1.set_names('x', 'PFC regional dimensions', ...
    'y', 'PFC predictive dimensions', ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g1.axe_property('XTickLabelRotation',35)
g1.draw()

%%
figure(233)
g2 = gramm(  'subset', subset,...
            'x', x2,...
            'y', y2,...
            'color', categorical(hpcsubset.patternAbstract),...
            'lightness', categorical(hpcsubset.control));
        
g2.facet_grid(categorical(hpcsubset.control), categorical(hpcsubset.patternAbstract))
g2.geom_point('dodge', 0.5, 'alpha', 0.3);
g2.geom_abline('style','k:');
g2.stat_cornerhist('edges',-4:0.5:4, 'aspect',1.2);
g2.set_point_options('base_size', 10);
g2.set_text_options('label_scaling', 1.5, 'base_size', 14);

g2.set_names('x', 'HPC regional dimensions', ...
    'y', 'HPC predictive dimensions', ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g2.axe_property('XTickLabelRotation',35)        
g2.draw()
%%
clf
figure(9)
g=gramm('x',cars.Horsepower,'y',cars.MPG,'subset',cars.Cylinders~=3 & cars.Cylinders~=5,'lightness',cars.Cylinders);

g(1,3).geom_point();

g(1,3).set_names('x','Horsepower','y','MPG','lightness','# Cyl');

g(1,3).set_title('lightness');
