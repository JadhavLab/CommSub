Patterns = Patterns_AllAnimals;

tmp = nd.fieldGet(Patterns,'rankRegress');

optDim = nd.fieldGet(tmp,'optDimReducedRankRegress');
optDim = ceil(median(optDim, 1));
optDim = max(optDim,[],'all');

% Run all of them
rt = table(); % Table for storing results


for iPatternType = progress(1:size(Patterns,1), 'Title','iPatternType')
for iPartition = progress(1:size(Patterns,2), 'Title','iPartition')
    for baseDirection = 1:size(Patterns,3)
    for removeDirection = 1:size(Patterns,3)
        for basePattern = 1:numel(patternnames)
        for removePattern = 1:numel(patternnames)
            
            % Base patterns
            generateH = Patterns(iPatternType, iPartition, baseDirection, basePattern).generateH;
            X_source  = Patterns(iPatternType, iPartition, baseDirection, basePattern).X_source;
            X_target  = Patterns(iPatternType, iPartition, baseDirection, basePattern).X_target;
            cvLoss    = Patterns(iPatternType, iPartition, baseDirection, basePattern).rankRegress.cvLoss;
            numDimsUsedForPrediction = 1:size(cvLoss,2);
            % Remove patterns
            B_ = Patterns(iPatternType, iPartition, removeDirection, removePattern).rankRegress.B_;
            % Run dim removal
            [performance,fullmodel] = plots.sequentialRemovePredDims(X_source, X_target, B_, optDim,...
                cvLoss, numDimsUsedForPrediction, "normalized", false);
            
            
            T = table(generateH, iPartition, baseDirection, removeDirection, patternnames(basePattern), patternnames(removePattern),fullmodel,...
                        'VariableNames',{'genH', 'partition','baseDirection','removeDirection','basePattern','removePattern','fullmodel'});
            T = repmat(T,numel(performance),1);
            T.dimensionRemoved = (0:numel(performance(:))-1)';
            T.performance = performance(:);
            T.performanceRemoved = T.performance(1) - T.performance;
            rt = [rt; T];
        end
        end
    end
    end
end
end

% Post processing of table data labels
C = @categorical;
rt.sameDirection = rt.baseDirection == rt.removeDirection;
dirlabel = ["remove other area","remove same area"]';
rt.sameDirectionLabel = dirlabel(rt.sameDirection+1);
dirlabel = ["hpc-hpc","hpc-pfc"]';
rt.targetArea = dirlabel(rt.baseDirection);
rt.removeTargetArea = dirlabel(rt.removeDirection);
rt.basePatternLabel = rt.basePattern.replace('theta','$\theta$').replace('delta','$\delta$').replace('ripple','SPW-R');
rt.removePatternLabel = rt.removePattern.replace('theta','$\theta$').replace('delta','$\delta$').replace('ripple','SPW-R');
rt.genHshort = shortcut.generateH(rt.genH);
