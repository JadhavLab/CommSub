function FigDat = plotexplained(FigDat, Option, varargin)
% plotPredictionPerPattern plots the prediction per pattern for the
% different methods and compares them with a KS-test to see if they are
% significantly different.

ip = inputParser;
ip.addParameter('figAppend', '', @(x) ischar(x) || isstring(x));
%ISSUE: maybe this no longer good range after zscore?
% ip.addParameter('xlim', [0,0.4]); % range or 'auto'
ip.addParameter('xlim', 'auto'); % range or 'auto'
ip.addParameter('nbins', 20); % number of bins or 'auto'
ip.addParameter('yscale', 'linear'); % 'linear' or 'log'
ip.parse(varargin{:});
Opt = ip.Results;

Option = munge.reshapeByLabels(Option, 1, [Option.genH_name]);

fig("predicition per pattern");
clf

[nMethods, nPatterns, ~] = size(FigDat.r_withhpc_patterns);

h_pred       = zeros(nMethods,nPatterns);
p_pred       = zeros(nMethods,nPatterns);
meanpred_hpc = zeros(nMethods,nPatterns);
meanpred_pfc = zeros(nMethods,nPatterns);

Nh = nPatterns;
Nw = nMethods;
gap = [0.01 0.03];
marg_h = [0.1 0.01];
marg_w = [0.01 0.01];
ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w);

for m = 1:nMethods
    for i = 1:nPatterns
        loc = (i-1)*nMethods + m;
        ax = ha(loc);
        axes(ax);

        curr1 = [FigDat.r_withhpc_patterns{m,i,:}];
        indices1 = intersect(find(curr1>0), find(~isnan(curr1)));
        hpcPreds = histogram(curr1(indices1), 'NumBins', Opt.nbins);
        set(hpcPreds, 'EdgeColor', 'none', 'FaceAlpha', 0.33);
        hold on
        
        curr2 = [FigDat.r_withpfc_patterns{m, i, :}];
        indices2 = intersect(find(curr2>0), find(~isnan(curr2)));
        pfcPreds = histogram(curr2(indices2), 'NumBins', Opt.nbins);
        set(pfcPreds, 'EdgeColor', 'none', 'FaceAlpha', 0.33);
        
        pattern_mean_withhpc = mean(FigDat.patternVarExplained_hpc(m,i,:));
        pattern_mean_withpfc = mean(FigDat.patternVarExplained_pfc(m,i,:));
        meanpred_hpc(m,i) = pattern_mean_withhpc;
        meanpred_pfc(m,i) = pattern_mean_withpfc;
        hold on

        ym = ylim();
        ym = ym(2);
        
        if m == 1
            avg_hh=line([pattern_mean_withhpc,pattern_mean_withhpc],[0 ym]);
            avg_hh.LineStyle = ':'; % Make line dotted
            avg_hh.LineWidth = 2;  % Thicken the line
            avg_hh.Color = 'blue';
            
            hold on
            avg_hp=line([pattern_mean_withpfc,pattern_mean_withpfc],[0 ym]);
            avg_hp.LineStyle = ':'; % Make line dotted
            avg_hp.LineWidth = 2;  % Thicken the line
            avg_hp.Color = 'red'; % Color it black
            
        else
            avg_hh=line([pattern_mean_withhpc,pattern_mean_withhpc],[0 ym]);
            avg_hh.LineStyle = ':'; % Make line dotted
            avg_hh.LineWidth = 2;  % Thicken the line
            avg_hh.Color = 'blue';
            
            hold on
            avg_hp=line([pattern_mean_withpfc,pattern_mean_withpfc],[0 ym]);
            avg_hp.LineStyle = ':'; % Make line dotted
            avg_hp.LineWidth = 2;  % Thicken the line
            avg_hp.Color = 'red'; % Color it black
        end
        
        
        if Option(1).sourceArea == "CA1"
            legend("hpc-hpc","hpc-pfc")
        else
            legend("pfc-hpc","pfc-pfc")
        end
        
        
        [h_pred(m,i),p_pred(m,i)] = kstest2(curr1(indices1), curr2(indices2));
        ylabel("data sets")
        xlabel("performance")
        pattern = FigDat.prophpc(m,i).name;
        genH = shortcut.generateH(FigDat.prophpc(m,i).generateH);
        title(pattern + " " + genH);
        if Opt.xlim ~= "auto"
            xlim(Opt.xlim) 
        end
    end
end

axs=findobj(gcf,'Type','axes');
if Opt.xlim == "auto"
    xlims=get(axs,'XLim');
    xlims=cat(1,xlims{:});
    xlims = [min(xlims(:,1)),max(xlims(:,2))];
    set(axs,'XLim',xlims)
end
set(axs, 'YScale', Opt.yscale)

folder = fullfile(codedefine, 'figures', 'new', '2b-prediction');
if ~exist(folder, 'dir')
    mkdir(folder)
end

set(gcf, 'Position', get(0, 'Screensize'));
savefig(strcat("Fig2b, prediction", Opt.figAppend, ".fig"))
saveas(gcf, fullfile(folder, strcat("Fig2b, prediction", Opt.figAppend, ".png")))
saveas(gcf, fullfile(folder, strcat("Fig2b, prediction", Opt.figAppend, ".pdf")))

FigDat.h_pred = h_pred;
FigDat.p_pred = p_pred;
FigDat.meanpred_hpc = meanpred_hpc;
FigDat.meanpred_pfc = meanpred_pfc;

end
