function plotStatOverX(b, hpc, pfc, stat, over, varargin)

ip = inputParser;
ip.addParameter('savedir', []);
ip.addParameter('patternSym',[]);
ip.addParameter('patternNames',[]);
ip.addParameter('qlim',[0.0, 0.95]);
ip.addParameter('plot_high_patterns_only', true);
ip.addParameter('normalize', false);
ip.addParameter('animal', []);
ip.parse(varargin{:})
Opt = ip.Results;
savedir = Opt.savedir;

if ~isempty(Opt.animal)
    for field = string(fieldnames(hpc.org))'
        hpc.org.(field) = hpc.org.(field)(a,:,:,:,:,:,:,:);
        pfc.org.(field) = pfc.org.(field)(a,:,:,:,:,:,:,:);
    end
end

[animals, partitions, phases, patterns, components, epochs] = size(hpc.org.sv);
if Opt.plot_high_patterns_only
    patterns = patterns/2;
end

switch char(stat)
    case "sv"
        statname = "variance";
    case "m"
        statname = "mean";
    case "am"
        statname = "abs mean";
    otherwise
        error();
end

titles = ["reward","error","outbound","inbound"];
directs = ["hpchpc", "pfchpc"];
screensize = get(0,'screensize');
screensize(3) = screensize(3)/3;
screensize(4) = screensize(4)/3;
alphaVal = 8/(animals*partitions);
for i = 1:2
    for j = 1:4
        f(i,j) = sfigure('Position',screensize);
        clf
        tiledlayout(patterns,components,'TileSpacing', 'compact');
        sfigure(f(i,j));
        titlestr = over + " " + directs(i) + " " + titles(j) + " : component " + statname;
        sgtitle(titlestr);
    end
end
cnt = 0;
totalSubplots = phases * patterns * components;
colors1 = crameri('buda', patterns * components);
colors2 = crameri('imola', patterns * components);
norm = @(x) (x - min(x))./(max(x) - min(x));
for a = progress(1:animals,'Title','paritions')
    for p = progress(1:partitions,'Title','paritions')
        for i = 1:phases
            for j = 1:patterns
                for k = 1:components
                    cnt = cnt + 1;
                    ind = sub2ind([components, patterns], k, j);
                    sfigure(f(1, i));
                    nexttile(ind);
                    hold on
                    Y = squeeze(hpc.org.(stat)(a, p, i, j, k,:));
                    if Opt.normalize
                        Y = norm(Y);
                    end
                    X = 1:numel(Y);
                    if over == "performance"
                        X = b.tperf(a, i, X);
                    end
                    s = scatter(X(:),Y(:));
                    set(s, 'MarkerFaceColor', colors1(ind,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alphaVal);
                    set(0, 'DefaultFigureRenderer', 'painters');
                    % print(gcf,'-depsc','-painters','myEps')
                    if p == 1
                        title(Opt.patternSym(j) + ", K\rightarrow" + num2str(k));
                        set(gca, 'color', [0.5,0.5,0.5]);
                        if j == patterns
                            xlabel(over)
                        end
                    elseif p == partitions
                        ydata = get(findobj(gca,'type','scatter'),'ydata');
%                         ylim = quantile([ydata{:}], Opt.qlim);
                         ylim = [0, 10];
                        set(gca,'ylim',ylim);
                    end
                    sfigure(f(2, i));
                    nexttile(ind);
                    Y = squeeze(pfc.org.(stat)(a, p, i, j, k,:));
                    X = 1:numel(Y);
                    if over == "performance"
                        X = b.tperf(a, i, X);
                    end
                    hold on;
                    if Opt.normalize
                        Y = norm(Y(:));
                    end
                    s = scatter(X(:),Y(:));
                    set(s, 'MarkerFaceColor', colors2(ind,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alphaVal);
                    set(0, 'DefaultFigureRenderer', 'painters');
                    %r = polyfit(X(:),Y(:),1);
                    if p == 1
                        title(Opt.patternSym(j) + ", K\rightarrow" + num2str(k));
                        set(gca, 'color', [0.5,0.5,0.5]);
                        if j == patterns
                            xlabel(over)
                        end
                    elseif p == partitions
                        ydata = get(findobj(gca,'type','scatter'),'ydata');
                      %  ylim = quantile([ydata{:}], Opt.qlim);
                      ylim = [0, 10];
                        set(gca,'ylim',ylim);
                    end
                    
                end
            end
        end
    end
end

% Median of measure
for i = progress(1:phases,'Title','Phases')
    for j = 1:patterns
        for k = 1:components
            cnt = cnt + 1;
            ind = sub2ind([components, patterns], k, j);
            sfigure(f(1, i));
            nexttile(ind);
            hold on
            if ~isempty(Opt.animal)
                Y = mean(squeeze(hpc.org.(stat)(:, :, i, j, k,:)), [1 2]);
            else
                Y = median(squeeze(hpc.org.(stat)(:, :, i, j, k,:)), [1 2]);
            end
            if Opt.normalize
                Y = norm(Y);
            end
            X = 1:numel(Y);
            if over == "performance"
                X = b.tperf(:, i, X);
            end
            s = plot(X(:),Y(:));
            set(s, 'Color', 'White');
            set(0, 'DefaultFigureRenderer', 'painters');
            %r = polyfit(X(:),Y(:),1);
            sfigure(f(2, i));
            nexttile(ind);
            if ~isempty(Opt.animal)
                Y = mean(squeeze(pfc.org.(stat)(:, :, i, j, k,:)),[1 2]);
            else
                Y = median(squeeze(pfc.org.(stat)(:, :, i, j, k,:)),[1 2]);
            end
            X = 1:numel(Y);
            if over == "performance"
                X = b.tperf(:, i, X);
            end
            hold on;
            if Opt.normalize
                Y = norm(Y(:));
            end
            s = plot(X(:),Y(:));
            set(0, 'DefaultFigureRenderer', 'painters');
            %r = polyfit(X(:),Y(:),1);
            set(s, 'Color', 'white');
        end
    end
end

for i = 1:2;
    for j = 1:4
        titlestr = over + " " + directs(i) + " " + titles(j) + " : component " + statname;
        if Opt.plot_high_patterns_only
            titlestr = titlestr + "_highOnly";
        end
        if Opt.normalize
            titlestr = titlestr + "_normalize";
        end
        export_fig titlestr.svg -painters;
        %saveas/export_fig(f(i,j), fullfile(savedir, titlestr + ".png"));
    end
end

