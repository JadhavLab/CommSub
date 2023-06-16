% run allPatternsSubspaceSimilarity.m first
% also need variables set up in component_Plots....
%% gramm plot of the relative average strength
fig('Mean strength gramm plot'); clf
patNames = replace(patternNames, 'control','low');
x = categorical(patNames);
y = rewardComponentPatternwise;

g = gramm('x', x ,'y', y);
g.geom_bar()
g.draw()

%% gramm plot of the relative average strength of a certain dimension
figure(112); clf
% 4*6 plot???
y = components.organizeComponents(rewardedComponent',5,3,"dim");

new_x = [];
for i = 1:size(y,1)
    new_x = [new_x;x];
end

figure;
g = gramm('x', new_x(:) ,'y', y(:));
g.geom_point()
g.stat_summary('type', 'std','geom','bar')
g.draw()

%%

