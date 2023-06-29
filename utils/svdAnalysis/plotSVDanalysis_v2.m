% Analyses for the v2 version of the singular value analysis
resFields = ["rotation","magnitude","magAdjRotation"];
measure = "rotation";
mfun = @(x) max(abs(quantile(x(:), [0.005, 0.995])));
qfun = @(x) [-mfun(x), mfun(x)];
analyses = struct();

% --------------------------------------------
% 1. R_n(p, c) : neuron wise pattern component
% --------------------------------------------

head(tab.U)
[groups, neuron, pattern, component] = findgroups(tab.U.neuron, tab.U.pattern, tab.U.component);
R = table(neuron, pattern, component);
for field = resFields
    R.(field) = splitapply(@mean, tab.U.(field), groups);
end

% ............................................
% Visualizing the pattern, component landscape
% ............................................
wrap = 8;
nNeuron = max(neuron);
fig(measure);
uNeuron = unique(neuron);
for N = progress(uNeuron')
    mat = unstack(R(R.neuron==N,["pattern","component",measure]), measure, "component", "VariableNamingRule","preserve");
    mat = table2array(mat(:,2:end));
    subplot(ceil(nNeuron/wrap), wrap, N); imagesc(mat);
    title("Neuron " + N);
    cmocean('balance');
    axis square
    analyses.neuronalComponents(N,:,:) = mat;
end
set(findobj(gcf,'type','axes'),'clim',qfun(analyses.neuronalComponents))

% ................................................
% Visualizing the ABS pattern, component landscape
% ................................................
nNeuron = max(neuron);
fig(measure + " abs");
uNeuron = unique(neuron);
for N = progress(uNeuron')
    mat = unstack(R(R.neuron==N,["pattern","component",measure]), measure, "component", "VariableNamingRule","preserve");
    mat = table2array(mat(:,2:end));
    subplot(ceil(nNeuron/wrap), wrap, N); imagesc(abs(mat));
    title("Neuron " + N);
    cmocean('balance');
    set(gca,'clim',[-0.15,0.15])
    axis square
    analyses.neuronalComponentsAbs(N,:,:) = mat;
end


% ............................
% Activation in + / - patterns
% ............................
nNeuron = max(neuron);
fig(measure + " abs +/i");
uNeuron = unique(neuron);
axstack=[];
clear area
for N = progress(uNeuron')
    mat = unstack(R(R.neuron==N,["pattern","component",measure]), measure, "component", "VariableNamingRule","preserve");
    mat = table2array(mat(:,2:end));
    mat_plus = sum(abs(mat(1:3,:)),1);
    mat_minus = sum(abs(mat(4:6,:)),1);
    axstack=[axstack, subplot(ceil(nNeuron/wrap), wrap, N)]; 
    hold off;
    A=area(1:size(mat_plus,2), mat_plus-mat_minus);
    set(A,'FaceAlpha',0.5,'FaceColor','k')
    hold on
    scatter(1:size(mat_plus,2), mat_plus-mat_minus, 'ko');
    hline(0)
    legend('plus - minus')
    title("Neuron " + N);
    cmocean('balance');
    set(gca,'clim',[-0.15,0.15])
    axis square
    analyses.neuronalComponentsAbsDiff(N,:,:) = mat_plus-mat_minus;
end
set(findobj(gcf, 'type','axes'), 'clim', analyses.neuronalComponentsAbsDiff(N,:,:));
linkaxes(axstack,'xy');

% ........................................
% Neuron wide - Patterns have more content
% ........................................
clear g
Rm = R(R.pattern >= 4,["pattern","component",measure]);
Rp = R(R.pattern <= 3, ["pattern","component",measure]);
Rdiff  = Rp;
Rdiff.(measure) = abs(Rp.(measure)) - abs(Rm.(measure));
analyses.gramm_neuronalAbsDiff(N,:,:) = Rdiff;
f=fig(measure + " abs +/i -- all neuron");close(f);
f=fig(measure + " abs +/i -- all neuron");
g = gramm('x', Rdiff.component, 'y', Rdiff.(measure));
g.set_color_options('chroma',0,'lightness',0)
g.stat_summary('geom',{'lines','area','point'})
g.geom_hline('yintercept', 0, 'style','r--');
g.geom_hline('yintercept', mean(Rdiff(Rdiff.component>2,:).(measure)), 'style','k:');
%g.axe_property('yscale','log')
g.draw()

% ........................................
% Neuron wide - Patterns have more content
% ........................................
clear g
Rm = R(R.pattern >= 4, ["pattern","component",measure]);
Rp = R(R.pattern <= 3, ["pattern","component",measure]);
Rdiff  = Rp;
Rdiff.(measure) = abs(Rp.(measure)) - abs(Rm.(measure));
f=fig(measure + " abs +/i -- all neuron");close(f);
f=fig(measure + " abs +/i -- all neuron");
g = gramm('x', Rdiff.component, 'y', Rdiff.(measure));
g.facet_grid(Rdiff.pattern, []);
%g.set_color_options('chroma',0,'lightness',0)
g.stat_summary('geom',{'lines','area','point'})
g.geom_hline('yintercept', 0, 'style','r--');
g.geom_hline('yintercept', mean(Rdiff(Rdiff.component>2,:).(measure)), 'style','k:');
%g.axe_property('yscale','log')
g.draw()

% -----------------------------------------------------------
% 2 R_n(p, c) * R_n(p, c).T : neuron wise pattern correlation
% -----------------------------------------------------------
wrap = 8;
nNeuron = max(neuron);
f=fig(measure + " pattern component correlation");
uNeuron = unique(neuron);
for N = progress(uNeuron')
    mat = R(R.neuron==N, ["pattern","component",measure]);
    mat = unstack(mat, measure, "component", "VariableNamingRule", "preserve");
    mat = table2array(mat(:,2:end));
    mat = mat * mat';
    for k = 1:size(mat,1); mat(k,k) = 0; end
    subplot(ceil(nNeuron/wrap), wrap, N); 
    imagesc(mat);
    title("Neuron " + N);
    cmocean('balance');
    axis square
    analyses.patternCorrelation(N,:,:) = mat;
end
set(findobj(gcf, 'type','axes'), 'clim', qfun(analyses.patternCorrelation(:,:,:)));

% --------------------------------------------------------------------------
% 3 Sum( R_n(p, c) * R_n(p, c).T ) : Sum of neuron-wise pattern correlations
% --------------------------------------------------------------------------
wrap = 8;
nNeuron = max(neuron);
f=fig(measure + " pattern component correlation average");
imagesc(squeeze(...
    median(abs(analyses.patternCorrelation(:,:,:)),1)...
    ));
cmocean('balance');
set(gca, 'clim', [-0.005,0.005])
%set(findobj(gcf, 'type','axes'), 'clim', qfun(analyses.patternCorrelation(:,:,:)));
