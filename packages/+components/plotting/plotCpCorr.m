
% reorganize the matrix in different ways and then plot imagesc of the
% correlation

for i = 1:4
    components_by_strength{i} = components.organizeComponents...
        (Patterns(1,2,1).compAnalysisByBeh(i).critical_components', 5, 3, "strength");
    components_by_dim{i} = components.organizeComponents...
        (Patterns(1,2,1).compAnalysisByBeh(i).critical_components', 5, 3, "dim");
    components_by_strength{i}(isnan(components_by_strength{i})) = 0;
    components_by_dim{i}(isnan(components_by_dim{i})) = 0;
    
end

% diff axis ticks!!
% by strength plot
fig('Correlation of Components during Critical Behaviors - by Pattern Strength'); clf; tiledlayout(2,2);
titles = ["reward", "error", "inBoundChoice","outBoundChoice"];
for i = 1:4
    nexttile;
   
    
    imagesc(corrcoef(components_by_strength{i}'))
    %cmocean haline
    crameri bam
    clim = get(gca,'clim');
    
    title(titles(i));
    yticks([0,5,10,15,20,25])
    yticklabels(["theta", "low-theta","delta","low-delta","ripple","low-ripple"])
end


% by dim plot
fig('Correlation of Components during Critical Behaviors - by # Dimension'); clf; tiledlayout(2,2);
titles = ["reward", "error", "inBoundChoice","outBoundChoice"];
for i = 1:4
    nexttile;
     R_bybeh{i} = corrcoef(components_by_dim{i}');
    imagesc(corrcoef(components_by_dim{i}'))
    %cmocean haline
    crameri bam
    clim = get(gca,'clim');
    
    title(titles(i));
    yticks([1,7,13,19,25])
    yticklabels(["dim 1", "dim 2","dim 3","dim 4","dim 5"])
    xticks([1,7,13,19,25])
    xticklabels(["dim 1", "dim 2","dim 3","dim 4","dim 5"])
end
%% All components' correlation
all_components_by_dim = components.organizeComponents(behavior_running_subspaces,5,3,"dim");
all_components_by_strength = components.organizeComponents(behavior_running_subspaces,5,3,"strength");

figure;
R = corrcoef(all_components_by_dim');
imagesc(corrcoef(all_components_by_dim')),crameri bam
yticks([3,9,15,21,27])
yticklabels(["dim 1", "dim 2","dim 3","dim 4","dim 5"])
xticks([3,9,15,21,27])
xticklabels(["dim 1", "dim 2","dim 3","dim 4","dim 5"])
axis square
    
    for j = 1:4
        lineObject = line([0 30] ,[j*6+0.5 j*6+0.5]);
        lineObject.LineStyle = ':'; % Make line dotted
        lineObject.LineWidth = 2;  % Thicken the line
        lineObject.Color = 'white'; %
        hold on
    end
        for j = 1:4
        lineObject = line([j*6+0.5 j*6+0.5] ,[0 30]);
        lineObject.LineStyle = ':'; % Make line dotted
        lineObject.LineWidth = 2;  % Thicken the line
        lineObject.Color = 'white'; %
        hold on
    end
colorbar

%% only the strong pattern actiity
all_cp_pattern = [];
for i = 1:3
    all_cp_pattern = [all_cp_pattern; all_components_by_dim((i-1)*6+1:(i-1)*6+3,:)];
end

figure;
R = corrcoef(all_cp_pattern');
imagesc(corrcoef(all_cp_pattern')),crameri bam
yticks([2,5,8,11,14])
yticklabels(["dim 1", "dim 2","dim 3","dim 4","dim 5"])
xticks([2,5,8,11,14])
xticklabels(["dim 1", "dim 2","dim 3","dim 4","dim 5"])
axis square
    
    for j = 1:4
        lineObject = line([0 30] ,[j*3+0.5 j*3+0.5]);
        lineObject.LineStyle = ':'; % Make line dotted
        lineObject.LineWidth = 3;  % Thicken the line
        lineObject.Color = 'white'; %
        hold on
    end
        for j = 1:4
        lineObject = line([j*3+0.5 j*3+0.5] ,[0 30]);
        lineObject.LineStyle = ':'; % Make line dotted
        lineObject.LineWidth = 3;  % Thicken the line
        lineObject.Color = 'white'; %
        hold on
    end
colorbar