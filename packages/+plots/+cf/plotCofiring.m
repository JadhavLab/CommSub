function FigDat = plot_cofiring(FigDat, Option, varargin)
% should use a two-sample t-test, since the co-firing of hpc and pfc are
% independently measured

ip = inputParser;
ip.addParameter('figAppend', '', @(x) ischar(x) || isstring(x));
ip.addParameter('normalization', 'count', @(x) ischar(x) || isstring(x));
ip.addParameter('quantileMarkers', [0.5], @(x) isnumeric(x) && isvector(x));
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
    plotline(FigDat, i, cofiring_hh, Opt, 'blue', "hpc-hpc corfiring mean")

    hold on
    plotline(FigDat, i, cofiring_hp, Opt, 'red', "hpc-pfc corfiring mean")
    
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

half_screenwidth = get(0,'ScreenSize');
half_screenwidth(3) = half_screenwidth(3)/2;
set(gcf, 'Position', half_screenwidth);
        
fullfigname = fullfile(codedefine,'figures', 'new','cofire',...
    string(Opt.normalization) + "cofiring per pattern" + Opt.figAppend);
saveas(gcf, fullfigname + ".png")
saveas(gcf, fullfigname + ".pdf")
savefig(gcf, fullfigname + ".fig")

function plotline(FigDat, i, cofiring_hh, Opt, color, displayname)
    if isequal(Opt.quantileMarkers, 'mean')
        x= [FigDat.mean_withhpccorr_pattern(i),...
            FigDat.mean_withhpccorr_pattern(i)];
        y= [0 max(cofiring_hh.Values)];
    elseif isnumeric(Opt.quantileMarkers)
        x= zeros(numel(Opt.quantileMarkers),2);
        for j = 1:numel(Opt.quantileMarkers)
            x(j,:)= quantile(FigDat.withhpccorr_pairs{i}, Opt.quantileMarkers(j));
        end
        x= [quantile(FigDat.withhpccorr_pairs{i}, Opt.quantileMarkers),...
            quantile(FigDat.withhpccorr_pairs{i}, Opt.quantileMarkers)];
        y= [0 max(cofiring_hh.Values)];
    else
        error("quantileMarkers must be 'mean' or numeric vector")
    end
    for i = 1:size(x,1)
        avg_hh(i)=line(x(i,:),y);
        avg_hh(i).LineStyle = ':'; % Make line dotted
        avg_hh(i).LineWidth = 2;  % Thicken the line
        avg_hh(i).Color = color;
        avg.hh(i).DisplayName = displayname;
    end
end
end
