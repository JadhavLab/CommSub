load('~/Downloads/multiGenH.mat');

Patterns = cellfun(@(x) shiftdim(x,-1), Patterns, 'UniformOutput', false);
Patterns = cat(1, Patterns{:});
Patterns = nd.dimLabel(Patterns, 1, "generateH",      string(Option.generateH)); % Label the generateH dimension of struct
Patterns = nd.dimLabel(Patterns, 2, "iPartition");                               % Label the generateH dimension of struct
Patterns = nd.dimLabel(Patterns, 3, "directionality", directionality);           % Label the directionality dimension of struct
Patterns = nd.dimLabel(Patterns, 4, "name",           patternNames);             % Label the directionality dimension of struct


%%---------------------------------------------------
%% 0 . Acquire data from pattern struct in workspace
%% --------------------------------------------------

T = query.getPatternTable(Patterns, Option);

patternType = string(T.patternType);
patternType( ~contains(patternType, '-' ), : ) = patternType( ~contains(patternType, '-' ), : ) + "-";
columns = patternType.split('-');
columns(columns(:,2)=="",2) = "pattern activity";
T.patternAbstract = columns(:,1);
T.control = columns(:,2);
T.patternAbstractSymbol = T.patternAbstract.replace('theta','$\theta$').replace('delta','$\delta$').replace('ripple','SPW-R');
T.genH = T.generateH.replace('  ','').replace('fromRipTimes','').replace('from','').replace('Wpli','WPLI').replace('FilteredEEG','Hilbert');
T.genDir = T.genH + " " + T.directionality;
% T.winSize = cell2mat(T.winSize);

highlow = true;
if highlow  % if control versus noncontrol is low versus high level of a pattern
    T.highlow = categorical(string(T.control).replace('control','low').replace("pattern activity","high"));
end

% ------------------------------
% Throw out columns that are nan
% ------------------------------
T(:,any(ismissing(T),1)) = [];


% Convert any string columns to categorical
T = util.table.string2categorical(T);

%%------------------------------------------------------------
%% I . Wpli versus Hilbert : Same partitions
%% ------------------------------------------------------------

% Options
field = 'rrDim';
corner_kws = {'edges',[-0.8:0.05:0.8]*6, 'aspect',1, 'location', [8], 'fill',...
'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM
%ARE SET MANUALLY! you need this to surround your data
% LOCATION SETS WHERE THE CORNER HISTS ARE PLACED
T = util.table.categorical2string(T);
% Select RR-dimension of hpc-pfc hilbert
filt = ["$directionality=='hpc-pfc'",...
        "contains($generateH, 'EEG')"];
[x, t] = query.query(T,filt,field);
% Select RR-dimension of hpc-pfc Wpli
filt = ["$directionality=='hpc-pfc'",...
        "contains($generateH, 'Wpli')"];
y = query.query(T,filt,field);

f = figure;
set(f, 'Position', get(0,'ScreenSize'));
g = gramm('x', x,...
          'y', y,...
          'subset', t.control == "control")
% assert(all(t.patternType == pfcsubset.patternAbstract))
g.facet_grid(categorical(t.patternAbstractSymbol), []);
g.geom_jitter('alpha', 0.01,...
              'width', 0.02,...
              'height',0.02);
g.stat_cornerhist(corner_kws{:});
g.set_point_options('base_size', 10);
g.set_text_options('label_scaling', 1.5, 'base_size', 10);
g.set_names('x', "HPC-PFC" +newline+ "Hilbert", ...
    'y', "HPC-PFC" +newline+ "WPLI", ...
    'Color', 'Pattern (Grey is low)', ...
    'row','',...
    'Lightness', 'Treatment/Control');

g = g.set_color_options('chroma',0);
g.set_text_options('interpreter','latex','base_size',10)
%set(g.results.geom_jitter_handle,'MarkerSize',5)
g.draw()
results1 = g.results;
g.update('subset', t.control ~= "control",...
         'color', categorical(t.patternAbstract));
g.set_color_options();
g.stat_cornerhist(corner_kws{:});
g.geom_abline('style','k:');
% g.axe_property('XTickLabelRotation',35, 'axis','square')
g.geom_jitter('alpha', 0.50, 'width',0.35, 'height',0.35);
g.draw()
%set(g.results.geom_jitter_handle,'MarkerSize',5);

% ................................
% Add median to each corner_hist
% ................................
[G, patAb, ctrl] = findgroups(t.patternAbstract, t.control);
[~,~,uPatAb] = unique(patAb);
[~,~,uCtrl] = unique(ctrl);
med_hpc = splitapply(@median, t.(field), G); % get median of each group
[G, patAb, ctrl] = findgroups(t.patternAbstract, t.control);
[~,~,uPatAb] = unique(patAb);
[~,~,uCtrl] = unique(ctrl);
med = splitapply(@median, x - y, G); % get median of each group
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

% .........................................
% Permuation tests : WPLI high versus low
% .........................................
[G, patAb] = findgroups(t.patternAbstract);
p = zeros(size(unique(G)));
for gg = unique(G)'
    t_pat = x(G == gg & t.control == "pattern activity") - y(G==gg & t.control=="pattern activity");
    t_ctrl = x(G == gg & t.control == "control") - y(G==gg & t.control=="control");
    p(gg) = permutationTest(t_pat, t_ctrl, 1000, 'plotresult', 1, 'meanfunc', @nanmedian);
    title(patAb(gg));
end


% ........................................
% Permuation tests : Hilbert versus wpli
% ........................................
[G, patAb] = findgroups(t.patternAbstract);
p = zeros(size(unique(G)));
for gg = unique(G)'
    t_pat = x(G == gg & t.control == "pattern activity") - y(G==gg & t.control=="pattern activity");
    t_ctrl = x(G == gg & t.control == "control") - y(G==gg & t.control=="control");
    p(gg) = permutationTest(t_pat, t_ctrl, 1000, 'plotresult', 1, 'meanfunc', @nanmedian);
    title(patAb(gg));
end

%%--------------------------------------------
%% II . Wpli versus Hilbert : Summary figure
%% -------------------------------------------
% ..............................
% Delta coherence : Average dimension bar graph
% ..............................
% First, we need to unstack the table in such a way such that measurements with
% different directionality and generateH are in the same row.
T = util.table.string2catetorical(T);
indvar  = ["directionality", "generateH"];
datavar = "rrDim";
g = gramm('x', T.genDir, 'y', T.rrDim, 'lightness', T.control, 'color', T.genH) %'subset', string(T.patternAbstract)=="delta");
g.facet_grid(T.patternAbstract, T.directionality)
g.stat_summary(...
    'geom',{'bar','black_errorbar'},...
    'width',3.5,...
    'dodge', 3.5)
figure;
g.draw()
set(findobj(g.facet_axes_handles,'type','patch'),'facealpha',0.5)
%[index, indicator_varaible] = query.getIndex(T, [indvar]);
%assert(~ismember(indvar, index));
%T.indicator_variable = indicator_variable;
%TT = T;
%%TT(:, [entangled_vars, query.getDataVars(T,datavar)]) = [];
%TT = unstack(TT, datavar, indvar(1), 'ConstantVariables', [index,query.getDataVars(T,datavar)]);
%
%% We want to sort boxes by both directionality and generateH
%TT.genDir = TT.


%%------------------------------------------------------------
%% II . Wpli high v. Wpli low v. Hilbert : Correlations
%% ------------------------------------------------------------

% We want a vector for every combination of directionality and generateH methods
T = util.table.categorical2string(T);
[groups, pattern, direct, genH, highlow] = findgroups(T.patternAbstract, T.directionality, T.genH, T.highlow);
uGroups = unique(groups);

vector = cell(numel(uGroups),1);
for g = uGroups'
    G = g == groups;

    % For each group, get a vector
    vector{g} = T.(field)(G);
end
vector = cat(2, vector{:});

%V = corrcoef(vector);
V = vector;
RowNames = join([pattern, direct,genH,highlow], "-").replace('Hilbert','Power');
correlation = array2table(V);
correlation.Properties.VariableNames = RowNames;
writetable(correlation, fullfile(datadefine, 'genH_correlations.csv'));

%Python has some nicer tools for correlation, so these are getting sent to an ipython notebook for analysis
