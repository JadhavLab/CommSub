function FigDat = plot_cofiring(FigDat, Option, varargin)
% should use a two-sample t-test, since the co-firing of hpc and pfc are
% independently measured

ip = inputParser;
ip.addParameter('figAppend', '', @(x) ischar(x) || isstring(x));
ip.addParameter('normalization', 'count', @(x) ischar(x) || isstring(x));
ip.addParameter('quantileMarkers', {0.5,0.75,"mean"}, @iscell);
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
    hh=plotline(FigDat, i, cofiring_hh, Opt, 'blue', "hpc-hpc corfiring mean", "hpc");

    hold on
    hp=plotline(FigDat, i, cofiring_hp, Opt, 'red', "hpc-pfc corfiring mean", "pfc");
    
    if Option(1).sourceArea =="CA1"
        legend([cofiring_hh;cofiring_hp; hh; hp],["hpc-hpc", "hpc-pfc", string(Opt.quantileMarkers), string(Opt.quantileMarkers)])
    else
        error("not implemented")
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
saveas(gcf,  fullfigname + ".png")
saveas(gcf,  fullfigname + ".pdf")
savefig(gcf, fullfigname + ".fig")

function avg_hh = plotline(FigDat, i, cofiring_hh, Opt, color, displayname, area)
    avg_hh = gobjects(numel(Opt.quantileMarkers),1);
    for q = 1:numel(Opt.quantileMarkers)
        if isequal(Opt.quantileMarkers{q}, "mean")
            x = [mean(FigDat.("with"+area+"_pairs"){i}),...
                 mean(FigDat.("with"+area+"_pairs"){i})];
            y = [0 max(cofiring_hh.Values)];
        else
            x = [quantile(FigDat.("with"+area+"_pairs"){i}, Opt.quantileMarkers{q}),...
                 quantile(FigDat.("with"+area+"_pairs"){i}, Opt.quantileMarkers{q})];
            y = [0 max(cofiring_hh.Values)];
        end
        avg_hh(q)=line(x,y);
        avg_hh(q).LineStyle = ':'; % Make line dotted
        avg_hh(q).LineWidth = 2;  % Thicken the line
        avg_hh(q).Color = color;
        avg_hh(q).DisplayName = displayname;
        if isequal(Opt.quantileMarkers{q}, "mean")
            % add a black circle to top of the line
            plot(x(1),y(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
        else
            alpha(avg_hh(q), Opt.quantileMarkers{q});
        end
    end
end
end
