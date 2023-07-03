% -------------
% Shortcut vars
% -------------
onlyCoh = true;
Patterns = Patterns_AllAnimals;
directionality = ["hpc-hpc","hpc-pfc"];
patternNames = ["theta", "delta", "ripples","theta-control", "delta-control", "ripples-control"];
genH           = shortcut.generateH(Option.generateH);
direct         = shortcut.directionality(directionality);
patternSymbols = shortcut.patternSymbols(patternNames, 2);
figureFolder   = fullfile(figuredefine, 'subspaceAngle');

% ----------------------------
% Yank out the relevant data
% ----------------------------
rankRegress = nd.fieldGet(Patterns,'rankRegress');
B_ = nd.fieldGetCell(rankRegress,'B_');
V  = nd.fieldGetCell(rankRegress,"V");
sourceIndex = nd.fieldGetCell(Patterns,'index_source');
targetIndex = nd.fieldGetCell(Patterns,'index_target');

% ----------------------------
% How many components?
% ----------------------------
K = 1:6;
for measurement = ["similarity" ,"dissimilarity"]
    for normalize = [false, true]
        for k = progress(K, 'Title','K')

            % ---------------------
            % Compute dim-reduced B
            % ---------------------
            indices = nd.indicesMatrixForm(rankRegress);
            B = cell(max(indices,[],1));
            for index = progress(indices','Title','Creating B_rrr')
                I = num2cell(index);
                thisB_ = B_{I{:}};
                thisV = V{ I{:} };
                B{I{:}} = thisB_(:,1:k) * thisV(:,1:k)';
            end


            % Munge the shapes used to compute distances here
            partitionDim = 2;
            genHDim = 1; dirDim = 3; patternDim = 4;
            partitions = unique(indices(:,partitionDim));
            indices = unique(indices(:,[genHDim,dirDim,patternDim]),'rows');
            indices = indices(indices(:,end)<4,:);
            
            if onlyCoh
                genHs = genH(indices(:, genHDim));
                genHfilter = find(genHs == "filtered");
                indices = indices(genHfilter, :);
            end
            
            subspaceDist = zeros(numel(partitions),size(indices,1),size(indices,1));
            rowVar = string([]);
            getInd = @(part,rest) [rest(1) part rest(2:end)];
            
            
            % -----------------------------
            % Compute angle of each subpace
            % -----------------------------
            % (REDO : Once per partition)
            for part = progress(partitions','Title','Partitions')
            for i = 1:size(indices,1)
                for j = 1:size(indices,1)
                    I1 = getInd(part,indices(i,:));
                    I2 = getInd(part,indices(j,:));
                    I1 = num2cell(I1);
                    I2 = num2cell(I2);
                    subspaceDist(part,i,j) = subspace(B{I1{:}}, B{I2{:}});
                end

                [x, ~, y, z] = deal(I1{:});
                index     = [direct(y), patternSymbols(z), genH(x)];
                rowVar(i) = join(string(index),'-');

            end
            end
            if measurement == "similarity"
                subspaceDist = squeeze(median(subspaceDist,1));
                subspaceDist = subspaceDist + abs(min(subspaceDist,[],'all'));
            elseif measurement == "dissimilarity"
                subspaceDist = squeeze(median(subspaceDist,1));
            else
                error("You spelled it wrong")
            end
            if normalize
                subspaceDist = (subspaceDist-min(subspaceDist,[],'all')) ./(max(subspaceDist,[],'all')-min(subspaceDist,[],'all'));
            else
                % NOTHING :)
            end
            rowVar = cellstr(rowVar);
            analysisName = sprintf('%s%ssubspaceDist_K=%d_norm=%d_measure=%s.mat',...
                figureFolder,  filesep, k, normalize, measurement);
            save(analysisName, 'subspaceDist', 'rowVar');

            % ----------------------------
            % Normalize and plot distances
            % ----------------------------
            f = figc("subpsace distances, " + analysisName);
            G = clustergram(subspaceDist, 'ColumnLabels', rowVar, 'RowLabels', rowVar, 'Colormap', cmocean('balance'));
%             set(G, 'Colormap', crameri('lajolla'))
%             set(G, 'Colormap', crameri('-grayC'))
            set(G, 'Colormap', crameri('acton'))
            sgtitle(f.Name)
            cgFig = findall(0,'Type','Figure','Tag','Clustergram'); %handle to clustergram figure
            cgFig(1).Children(end).YTickLabelRotation=35;
            cgFig(1).Children(end).XTickLabelRotation=-25;
            cgFig(1).Children(end).FontSize=16;
            saveas(cgFig(1),fullfile(figureFolder,sprintf('subspaceDistance_K=%d.svg',k)));

        end
    end
end

%ax = findobj(gcf,'type','axes');
%set(ax,'YTickLabelRotation',35, 'XTickLabelRotation', -35, 'FontSize', 20);
%saveas(gcf,fullfile(figureFolder,'subspaceDistance.svg'))
