% Load data from your favorite BasicPatternFind into memory
% and run this script

%% Acquire tables for plotting
%[behWtab, tabs, vars] = characterize.makeBehaviorTable(data);
if isfield(data, 'folder')
    folder = data.folder;
else
    folder = ''
end

%tab = behWtab.trajdist.d;
%tab=sortrows(stack(tab, setdiff(tab.Properties.VariableNames, {'T','K'}), 'IndexVariableName', 'trajdist', 'NewDataVariableName', 'mu'),["K","T","trajdist"]);
%T=splitapply(@(x,y,z,zz) removevars(table(mean(x),mean(y),mean(z),mean(zz), 'VariableNames', tab.Properties.VariableNames),'T'), tab, findgroups(tab.K, tab.trajdist));
T = tabs{2};
T = removevars(T, ["RawL","stderr","P", "T"]);
G = findgroups(T.K, T.PrettyL);
TT = splitapply(@characterize.tab_meanAndMode, T, findgroups(T.K, T.PrettyL));


%% Pattern-wise label trends
%set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
fig(['Pattern-wise trends: ', folder])
clf
g = gramm('x', T.RawL, 'y', T.mu, 'color', T.K)
g.facet_grid(T.K, [], 'scale','free_y');
g.set_text_options('Interpreter','latex')
g.set_names('row','K','column','','x','Trajdist','y','Mu','color','K');
g.geom_bar()
g.draw()
%set(findobj(gcf,'type','axes'),'TickLabelInterpreter','latex')

%% Pattern-wise label trends
%set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
fig(['Pattern-wise trends: ', folder])
clf
q = quantile(T.mu, 0.95);
g = gramm('x', T.RawL, 'y', T.d, 'color', T.K, 'subset', T.mu > q)
g.facet_grid(T.K, [], 'scale','free_y');
g.set_text_options('Interpreter','latex')
g.set_names('row','K','column','','x','Trajdist','y','Mu','color','K');
%g.stat_bin('nbins', numel(unique(T.PrettyL)), 'geom', 'bar', 'normalization','pdf')
%g.geom_jitter('alpha',0.1)
g.stat_summary('type','sem', 'geom', {'area', 'black_errorbar'})
g.stat_surrogate_pvals()
g.draw()
set(findobj(gcf,'type','axes'),'TickLabelInterpreter','latex')


%% Pattern-wise label trends
%set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
fig(['Pattern-wise trends: ', folder])
clf
g = gramm('x', T.P, 'y', T.mu, 'color', T.K)
g.facet_grid(T.K, [], 'scale','free_y');
g.set_text_options('Interpreter','latex')
g.set_names('row','K','column','','x','Trajdist','y','Mu','color','K');
%g.stat_bin('nbins', numel(unique(T.PrettyL)), 'geom', 'bar', 'normalization','pdf')
%g.geom_jitter('alpha',0.1)
g.geom_jitter()
%g.stat_surrogate_pvals()
g.draw()
set(findobj(gcf,'type','axes'),'TickLabelInterpreter','latex')

