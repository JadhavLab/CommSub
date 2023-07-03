function [h_pred,p_pred,meanpred_hpc,meanpred_pfc] = plotVarExplained(FigDat, Option, varargin)
% plotPredictionPerPattern plots the prediction per pattern for the
% different methods and compares them with a KS-test to see if they are
% significantly different.

ip = inputParser;
ip.addParameter('figAppend', '', @(x) ischar(x) || isstring(x));
ip.parse(varargin{:});
Opt = ip.Results;

fig("predicition per pattern");
clf

[nMethods, nPatterns, ~] = size(FigDat.r_withhpc_patterns);

h_pred       = zeros(nMethods,nPatterns);
p_pred       = zeros(nMethods,nPatterns);
meanpred_hpc = zeros(nMethods,nPatterns);
meanpred_pfc = zeros(nMethods,nPatterns);

for m = 1:1
    for i = 1:nPatterns
        subplot(3,1,i)
        curr1 = [FigDat.r_withhpc_patterns{m}{i}{:}];
        indices1 = intersect(find(curr1>0), find(~isnan(curr1)));
        hpcPreds = histogram(curr1(indices1));
        set(hpcPreds, 'EdgeColor', 'none', 'FaceAlpha', 0.33);
        hold on
        
        curr2 = [FigDat.r_withpfc_patterns{m}{i}{:}];
        indices2 = intersect(find(curr2>0), find(~isnan(curr2)));
        pfcPreds = histogram(curr2(indices2));
        set(pfcPreds, 'EdgeColor', 'none', 'FaceAlpha', 0.33);
        
        pattern_mean_withhpc = mean(FigDat.patternVarExplained_hpc(m,i,:));
        pattern_mean_withpfc = mean(FigDat.patternVarExplained_pfc(m,i,:));
        meanpred_hpc(m,i) = pattern_mean_withhpc;
        meanpred_pfc(m,i) = pattern_mean_withpfc;
        hold on
        
        if m == 1
            avg_hh=line([pattern_mean_withhpc,pattern_mean_withhpc],[0 2000]);
            avg_hh.LineStyle = ':'; % Make line dotted
            avg_hh.LineWidth = 2;  % Thicken the line
            avg_hh.Color = 'blue';
            
            hold on
            avg_hp=line([pattern_mean_withpfc,pattern_mean_withpfc],[0 2000]);
            avg_hp.LineStyle = ':'; % Make line dotted
            avg_hp.LineWidth = 2;  % Thicken the line
            avg_hp.Color = 'red'; % Color it black
            
        else
            avg_hh=line([pattern_mean_withhpc,pattern_mean_withhpc],[0 2000]);
            avg_hh.LineStyle = ':'; % Make line dotted
            avg_hh.LineWidth = 2;  % Thicken the line
            avg_hh.Color = 'green';
            
            hold on
            avg_hp=line([pattern_mean_withpfc,pattern_mean_withpfc],[0 2000]);
            avg_hp.LineStyle = ':'; % Make line dotted
            avg_hp.LineWidth = 2;  % Thicken the line
            avg_hp.Color = 'yellow'; % Color it black
        end
        
        
        if Option.sourceArea == "CA1"
            legend("hpc-hpc","hpc-pfc")
        else
            legend("pfc-hpc","pfc-pfc")
        end
        
        
        [h_pred(m,i),p_pred(m,i)] = kstest2(curr1(indices1), curr2(indices2));
        ylabel("data sets")
        xlabel("performance")
        title(Option(1).patternNames(i));
        xlim([0,0.4])
        
    end
end

folder = fullfile(codedefine, 'figures', 'new', '2b-prediction');
savefig(strcat("Fig2b, prediction", Opt.figAppend, ".fig"))
saveas(gcf, fullfile(folder, strcat("Fig2b, prediction", Opt.figAppend, ".png")))
saveas(gcf, fullfile(folder, strcat("Fig2b, prediction", Opt.figAppend, ".pdf")))

FigDat.h_pred = h_pred;
FigDat.p_pred = p_pred;
FigDat.meanpred_hpc = meanpred_hpc;
FigDat.meanpred_pfc = meanpred_pfc;

end
