%{
-------------------------------
Possible subpsace angle methods
-------------------------------

1. **Principal Angles (Canonical Angles or Subspace Angles):** This method is
probably the most common way to measure the angle between two subspaces.
Principal angles are invariant under unitary transformations. In the context of
linear algebra, the principal angles are defined as the angles between pairs of
nearest points of two subspaces. If we have two subspaces \( U \) and \( V \)
of \( \mathbb{R}^n \) (or \( \mathbb{C}^n \)) with dimensions \( p \) and \( q
\) respectively, there are \( min(p, q) \) principal angles between \( U \) and
\( V \). The cosine of these angles can be found as the singular values of the
cross-covariance matrix \( X^T Y \).

2. **Frobenius Norm of Difference of Projections:** Another way to define the
angle between two subspaces is to look at the difference between the projection
matrices for each subspace. The projection matrix of a subspace \( S \)
represented by a matrix \( X \) can be calculated as \( X (X^T X)^{-1} X^T \).
The angle can then be defined as the Frobenius norm of the difference of the
projection matrices for the two subspaces.

3. **Geodesic Distance on the Grassmann Manifold:** This approach considers the
space of all \( k \)-dimensional subspaces in \( \mathbb{R}^n \) (or \(
\mathbb{C}^n \)) as a manifold (the Grassmann manifold) and defines the angle
between two subspaces as the geodesic distance between the two points on the
manifold representing the subspaces. This distance can be calculated from the
principal angles.

Note that these methods may give different results because they define the
angle in slightly different ways. The best method to use depends on the
specific context and requirements of the problem you are trying to solve.

%}


% -------------
% Shortcut vars
% -------------
onlyCoh = true;
directionality = ["hpc-hpc","hpc-pfc"];
patternNames = ["theta", "delta", "ripples","theta-control",...
                "delta-control", "ripples-control"];
genH           = shortcut.generateH([Option.generateH]);
direct         = shortcut.directionality(directionality);
patternSymbols = shortcut.patternSymbols(patternNames, 2);
figureFolder   = fullfile(figuredefine, 'subspaceAngle');

% Reshape so that animals in partition dimension 
% (they function as replicates in this analysis)
[n.animGenH, n.partition, n.dir, n.sing, n.patterns] = size(Patterns);
n.animals = length(unique([Patterns.animal]));
n.genH = length(unique([Patterns.genH_name]));
P=reshape(Patterns, [n.animals,n.genH,n.partition,n.dir,n.patterns]);
assert(length(unique([P(1,:,:,:,:,:).animal])) == 1, "Reshape failed")
disp("Size of Patterns: " + join(string(size(Patterns)), "x"))

% ----------------------------
% Yank out the relevant data
% ----------------------------
rankRegress = nd.fieldGet(P,'rankRegress');
rankRegress = ndb.toNd(rankRegress);
B_ = nd.fieldGetCell(rankRegress,'B_');
V  = nd.fieldGetCell(rankRegress,"V");
sourceIndex = nd.fieldGetCell(P,'index_source');
targetIndex = nd.fieldGetCell(P,'index_target');

% ----------------------------
% How many components?
% ----------------------------
K = 1:6;
for measurement = ["princ_similarity" ,"princ_dissimilarity", "frobenius_similarity", "frobenius_dissimilarity"]
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
                thisV = V{ I{:} }; % V is the same for all indices
                B{I{:}} = thisB_(:,1:k) * thisV(:,1:k)'; % B = B_ * V' (dim-reduced)
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
            if startsWith(measurement,"princ")
                ang_method = 'principal';
            elseif startsWith(measurement,"frobenius")
                ang_method = 'frobenius';
            else
                error("You spelled it wrong")
            end
            for part = progress(partitions','Title','Partitions')
                for i = 1:size(indices,1)
                    for j = 1:size(indices,1)
                        I1 = getInd(part,indices(i,:));
                        I2 = getInd(part,indices(j,:));
                        I1 = num2cell(I1);
                        I2 = num2cell(I2);
                        subspaceDist(part,i,j) = subspace(B{I1{:}}, B{I2{:}}, ang_method);
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
                subspaceDist = (subspaceDist-min(subspaceDist,[],'all')) ...
                ./(max(subspaceDist,[],'all')-min(subspaceDist,[],'all'));
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
