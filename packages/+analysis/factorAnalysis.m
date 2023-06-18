function Patterns = factorAnalysis(Patterns, Option)

consts = option.constants();
HPC = consts.HPC;
PFC = consts.PFC;

cvNumFolds = Option.rankRegress.cvNumFolds;
cvOptions = statset('crossval');
regressMethod = @FactorRegress;
for p = 1:Option.numPartition
    for i = 1:Option.nPatternAndControl
        for j = [HPC, PFC]
            q = 1: nSource;
            disp("processing the "+p+" partition and the "+i+" pattern "+j)
            currSource = Patterns(p,j,i).X_source';
            currTarget = Patterns(p,j,i).X_target';
            warning('')
            
            cvLoss = CrossValFa(currSource, q, cvNumFolds, cvOptions);
            qOpt =  FactorAnalysisModelSelect(cvLoss, q);
            
            cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
                (regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
                numDimsUsedForPrediction, ...
                'LossMeasure', 'NSE', 'qOpt', qOpt);
            
            cvl = crossval(cvFun, ...
                currTarget, currSource, ...
                'KFold', cvNumFolds, ...
                'Options', cvOptions);
            
            cvLoss = ...
                [ mean(cvl); std(cvl)/sqrt(cvNumFolds) ];
            
            Patterns(p,j,i).factorAnalysis.optDimFactorRegress = ModelSelect...
                (cvLoss, numDimsUsedForPrediction);
            
            Patterns(p,j,i).factorAnalysis.qOpt = qOpt;
            Patterns(p,j,i).factorAnalysis.cvl = cvl;
            Patterns(p,j,i).factorAnalysis.cvLoss = cvLoss;
            
            
            [Patterns(p,j,i).factorAnalysis.Z, ...
                Patterns(p,j,i).factorAnalysis.U,...
                Patterns(p,j,i).factorAnalysis.Q] = ...
                ExtractFaLatents(currSource, qOpt);
            [warnMsg, warnId] = lastwarn();
            
            if isempty(warnId)
                Patterns(p,j,i).singularWarning = false;
                
            else
                Patterns(p,j,i).singularWarning = true;
                disp("singular is marked")
            end
        end
    end
end
end

%%

for p = 1:Option.numPartition
for i = 1:nPatternAndControl
    for j = [HPC, PFC]
        try
            [Patterns(p,j,i).rankRegress.removedPerformance, ~,~]=...
                plots.getUncorrelatedPerformance( Patterns(p,j,i).rankRegress.B_,...
                Patterns(p,j,i).X_source, Patterns(p,j,i).X_target, Patterns(p,j,i).rankRegress.optDimReducedRankRegress,...
                numDimsUsedForPrediction, Patterns(p,j,i).rankRegress.cvLoss);
        catch
            keyboard
        end
        [Patterns(p,j,i).perfAsRemovingEach] = plots.sequentialRemovePredDims(Patterns(p,j,i).X_source, ...
            Patterns(p,j,i).X_target,Patterns(p,j,i).rankRegress.B_, Patterns(p,j,i).rankRegress.optDimReducedRankRegress, ...
            Patterns(p,j,i).rankRegress.cvLoss,numDimsUsedForPrediction);
    end
end

