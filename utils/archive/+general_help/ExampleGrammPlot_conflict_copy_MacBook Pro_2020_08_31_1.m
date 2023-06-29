%%

patternType = T.patternType;
patternType( ~contains(patternType, '-' ), : ) = patternType( ~contains(patternType, '-' ), : ) + "-";
columns = patternType.split('-');
columns(columns(:,2)=="",2) = "pattern activity";
T.patternAbstract = columns(:,1);
T.control = columns(:,2);

%% Average pattern maxDimPerc, single animal, single direction, all patterns
clf
figure(4)
winsizes = cell2mat(T.winSize);
subset = (T.animal == "JS15" | T.animal =="JS21") & T.directionality == "hpc-hpc" & T.generateH == "fromCoherence  fromRipTimes" & winsizes(:,2) == 0.125;

hpcsubset = T.directionality == 'hpc-hpc';
pfcsubset = T.directionality == 'hpc-pfc';

% x = categorical(T(hpcsubset, "rrDim"));
% y = T(pfcsubset, "rrDim");
x = categorical(T.patternType);
% y = T.percMax_rrDim; 
y = T.full_model_performance;
% g = gramm(  'x', x,...
%     'y', y);
% g.geom_point();
% g.stat_cornerhist('edges',-4:0.2:4,'aspect',0.6);
g = gramm('subset', subset, ...
    'x', x ,...
    'y', y,...
    'Color', categorical(T.patternAbstract), ...
    'lightness', categorical(T.control));

g.geom_jitter('alpha', 0.5, 'width', 0.6)
g.stat_summary('type', 'sem','geom','black_errorbar')
%g.stat_summary('type', 'sem','geom','bar')
g.set_point_options('base_size', 20);
g.set_text_options('label_scaling', 1.5, 'base_size', 14);

g.set_names('x', 'Pattern Type', ...
    'y', '#(Predictive Dims)/#(Max. Pred. Dims)', ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.set_title(" Coherence Example Pattern Dimensionalitiy" + newline() + "JS15 + JS21, HPC-PFC");

g.axe_property('XTickLabelRotation',35)
g.draw()

%% Cornerhist hpc versus pfc predictive
clf
figure(4)
winsizes = cell2mat(T.winSize);
%subset = (T.animal == "JS15" | T.animal =="JS21") & T.directionality == "hpc-hpc" & T.generateH == "fromCoherence  fromRipTimes" & winsizes(:,2) == 0.125;

hpcsubset = T.directionality == 'hpc-hpc';
pfcsubset = T.directionality == 'hpc-pfc';
hpcsubset = T(hpcsubset, :);
pfcsubset = T(pfcsubset, :);

x = hpcsubset.rrDim;
y = T(pfcsubset, "rrDim");
g = gramm(  'x', x,...
            'y', y);
g.facet_wrap(categorical(T.patternType))
g.geom_point();
g.stat_cornerhist('edges',-4:0.2:4, 'aspect',0.6);
g.set_point_options('base_size', 20);
g.set_text_options('label_scaling', 1.5, 'base_size', 14);
%g.set_title();
g.axe_property('XTickLabelRotation',35)
g.draw()
