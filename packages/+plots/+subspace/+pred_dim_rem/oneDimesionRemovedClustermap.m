figureFolder   = fullfile(figuredefine, 'subspaceAngle');
K = 2;

% Select all where K dimension removed
% ------------------------------------
RT = rt(t.dimensionRemoved == K, :);

% Group dimensions
% ----------------
x_dims = ["baseDirection", "basePattern"];
y_dims = ["removeDirection", "removePattern"];
possible_index_dims = ["genH", "partition", "baseDirection", "removeDirection", "basePattern", "removePattern"];
index_dims = setdiff(possible_index_dims, "partition");
data_dims = ["performanceRemoved"];

% Mean over any properties?
% -------------------------
if ~any(ismember(possible_index_dims,index_dims))

    I = num2cell(index_dims);
    [groups, index_vals{:}] = findgroups(I{:});
    uGroups = unique(groups);
    for g = uGroups'
        ivals = cellfun(@(x) x(g), index_vals,'UniformOutput',false);
        iField = 0;
        for field = data_dims
            iField = iField + 1;
            dvals{iField} = mean(RT(groups == g,:).(field));
        end
        tab = table(ivals{:}, dvals{:}, ...
            'VariableNames', [index_dims, data_dims]);
        newRT = [newRT; tab];
    end

    RT = newRT;
    clear newRT;
end

% Let's unstack into a distance matrix based on performance as a distance metric
[subspaceDist, rowVar] = sequentialDistanceMatrix(RT, index_dims, x_dims, y_dims);

% ----------------------------
% Normalize and plot distances
% ----------------------------
f=figc("dimremove distances, K = " + K);
subspaceDist = (subspaceDist-min(subspaceDist,[],'all')) ./(max(subspaceDist,[],'all')-min(subspaceDist,[],'all'));
G = clustergram(subspaceDist, 'ColumnLabels', rowVar, 'RowLabels', rowVar, 'Colormap', cmocean('balance'));
sgtitle(f.Name)
cgFig = findall(0,'Type','Figure','Tag','Clustergram'); %handle to clustergram figure
cgFig(1).Children(end).YTickLabelRotation=35;
cgFig(1).Children(end).XTickLabelRotation=-25;
cgFig(1).Children(end).FontSize = 14;

saveas(cgFig,fullfile(figureFolder,'subspaceDistance_fromDimRemoval.svg'))
saveas(cgFig,fullfile(figureFolder,'subspaceDistance_fromDimRemoval.png'))
