numUsedForPrediction = min(nTarget,nSource);
% make the averaged version
curr_cvLoss_rr = cell(10,2,3);
curr_cvLoss_fa = cell(10,2,3);
curr_rrDim = cell(10,2,3);
curr_qOptDim = cell(10,2,3);

for p = 1:10
    for i = 1:nPatterns
        for j = 1:2
            if ~Patterns(p,j,i).singularWarning
                curr_cvLoss_fa{p,j,i} = Patterns(p,j,i).factorAnalysis.cvLoss;
                curr_qOptDim{p,j,i} = Patterns(p,j,i).factorAnalysis.optDimFactorRegress;
                curr_cvLoss_rr{p,j,i} = Patterns(p,j,i).rankRegress.cvLoss;
                curr_rrDim{p,j,i} = Patterns(p,j,i).rankRegress.optDimReducedRankRegress;
            end
        end
    end
end

figure(2333)
clf
full_model_performance_rr = [];
full_model_performance_fa = [];

for i = 1:nPatterns
    for j = 1:2
        subplot(3,2,2*(i-1)+j)
        
        full_model_rr = plots.plotPredictiveDimensions(numUsedForPrediction,...
            curr_cvLoss(:,j,i), "optDim", curr_rrDim(:,j,i),"mode", "rr");
        
        
        full_model_performance_rr = [full_model_performance_rr, full_model_rr];
        xlim([0,12.5])
        hold on
        
        plot(1, full_model_rr,'^');
        ylim([0, full_model_rr+0.1])
        
        
        
        hold on
        
        full_model_fa = plots.plotPredictiveDimensions(numUsedForPrediction,curr_cvLoss(:,j,i),'optDim',curr_qOptDim, ...
           "mode", "fa", "color", "blue");
 
        full_model_performance_fa = [full_model_performance_fa, full_model_fa];
        xlim([0,12.5])
        hold on
        
        plot(1, full_model_fa,'*');
        ylim([0, full_model_fa+0.1])
        if j == 1
            ax1 = gca;
        else
            ax2 = gca;
        end
        title([Patterns(p,j,i).name Patterns(p,j,i).directionality])
    end
    linkaxes([ax1,ax2],'y')
    
end