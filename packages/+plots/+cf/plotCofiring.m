function FigDat = plot_cofiring(FigDat, Option, varargin)
% should use a two-sample t-test, since the co-firing of hpc and pfc are
% independently measured

ip = inputParser;
ip.addParameter('figAppend', '', @(x) ischar(x) || isstring(x));
ip.addParameter('normalization', 'count', @(x) ischar(x) || isstring(x));
ip.parse(varargin{:});
Opt = ip.Results;

if ~startsWith(Opt.figAppend, ' ')
    Opt.figAppend = " " + Opt.figAppend;
end

nPatterns = numel(FigDat.withhpc_pairs);
patternnames = Option(1).patternNames;

fig("cofiring per pattern" + Opt.figAppend);

FigDat.h_corrdiff = zeros(1,nPatterns);
FigDat.p_corrdiff = zeros(1,nPatterns);

for i = 1:nPatterns
    subplot(nPatterns,1,i);
    hold off;
    cofiring_hh = histogram(FigDat.withhpc_pairs{i}, 'Normalization', Opt.normalization);
    set(cofiring_hh, 'EdgeColor', 'none', 'FaceAlpha', 0.33);
    hold on
    cofiring_hp = histogram(FigDat.withpfc_pairs{i}, 'Normalization', Opt.normalization);
    set(cofiring_hp, 'EdgeColor', 'none', 'FaceAlpha', 0.33);
    
    title(patternnames(i))
    ylabel("pairs")
    
    hold on
    x= [FigDat.mean_withhpccorr_pattern(i),FigDat.mean_withhpccorr_pattern(i)];
    y= [0 max(cofiring_hh.Values)];
    avg_hh=line(x,y);
    avg_hh.LineStyle = ':'; % Make line dotted
    avg_hh.LineWidth = 2;  % Thicken the line
    avg_hh.Color = 'blue';
    avg.hh.DisplayName = "hpc-hpc corfiring mean";
%     
    hold on
    x = [FigDat.mean_withpfccorr_pattern(i), ...
         FigDat.mean_withpfccorr_pattern(i)];
    y = [0 max(cofiring_hh.Values)];
    avg_hp=line(x,y);
    avg_hp.LineStyle = ':'; % Make line dotted
    avg_hp.LineWidth = 2;  % Thicken the line
    avg_hp.Color = 'red'; 
    avg.hh.DisplayName = "hpc-pfc corfiring mean";
    
    if Option(1).sourceArea =="CA1"
        legend("hpc-hpc","hpc-pfc")
    else
        legend("pfc-hpc","pfc-pfc")
    end
    
    xlabel("pairwise correlation")
    groupA = double(FigDat.withhpc_pairs{i}(~isnan(FigDat.withhpc_pairs{i})));
    groupB = double(FigDat.withpfc_pairs{i}(~isnan(FigDat.withpfc_pairs{i})));
    [FigDat.h_corrdiff(i),FigDat.p_corrdiff(i)] = ttest2(groupA, groupB);
    xlim([-0.1,0.4])
end
        
fullfigname = fullfile(codedefine,'figures', 'new','cofire',...
    string(Opt.normalization) + "cofiring per pattern" + Opt.figAppend);
saveas(gcf, fullfigname + ".png")
saveas(gcf, fullfigname + ".pdf")
savefig(gcf, fullfigname + ".fig")

