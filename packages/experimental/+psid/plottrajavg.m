function varOfD = plottrajavg(avg, varargin)

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('dim',           [1, 2, 3]);
ip.addParameter('crameri',       'berlin');
ip.addParameter('starttile',     true);
ip.addParameter('zscore',        false);
ip.addParameter('leftright',     "together");
ip.addParameter('sortDimVar',    false);
ip.addParameter('dim4asSize',    []); % if not is empty, then needs to be a size range
ip.addParameter('showDirection', true); % if not is empty, then needs to be a size range
ip.addParameter('distinguishRight', true); % if not is empty, then needs to be a size range
ip.addParameter('interpolate',   0); % if not is empty, then needs to be a size range
ip.addParameter('drawnow',       false); % if not is empty, then needs to be a size range
ip.parse(varargin{:})
Opt = ip.Results;
dim = Opt.dim;

if Opt.zscore
    groups = findgroups(avg.epoch, avg.col);
    uGroups = unique(groups);
    for group = uGroups'
        if sum(group==groups)
            avg(group == groups, :).val = zscore(avg(group == groups, :).val);
        end
    end
end

[mx, Mx] = deal(min(avg(avg.col==dim(1),:).val), max(avg(avg.col==dim(1),:).val));
[my, My] = deal(min(avg(avg.col==dim(2),:).val), max(avg(avg.col==dim(2),:).val));
[mz, Mz] = deal(min(avg(avg.col==dim(3),:).val), max(avg(avg.col==dim(3),:).val));


if Opt.interpolate == 0
    nColors = numel(unique(avg.ldistbin));
else
    nColors = Opt.interpolate;
end
colors = crameri(Opt.crameri, nColors);
[groups, epochs, trajbound, rewarded, leftright] = findgroups(avg.epoch, avg.trajbound, avg.rewarded, avg.leftright);
uGroups = unique(groups);
uEpochs = unique(epochs);
nEpochs = numel(uEpochs);
if Opt.leftright == "together"
    nCond = numel(unique(trajbound)) * numel(unique(rewarded));
elseif Opt.leftright == "separate"
    nCond = numel(unique(trajbound)) * numel(unique(rewarded)) * numel(unique(leftight));
else
    error("Bad input")
end

if Opt.starttile
    t = tiledlayout(nCond, nEpochs, 'TileSpacing', 'tight');
end

sortrew = @(x) abs(1-x);

varOfD = nan(max(epochs), max(trajbound)+1, max(rewarded)+1, max(leftright)+1, numel(unique(avg.col)));

for group = progress(uGroups','Title','Plotting')
    if Opt.leftright  == "together"
        i = sub2ind([nEpochs, nCond], epochs(group)/2, trajbound(group) * 2 + sortrew(rewarded(group))+1);
    elseif Opt.leftright == "separate"
        i = sub2ind([nEpochs, nCond], epochs(group)/2, trajbound(group) * 2^2 + sortrew(rewarded(group)) * 2 + leftright(group) + 1);
    end
    nexttile(i);
    D = avg(groups == group, ["ldistbin","col", "val"]);
    D = unstack(D, 'val', 'col', 'NewDataVariableName', "Dim" + (1:numel(unique(D.col))));

    D = table2array(D);
    if Opt.sortDimVar
        [~,idx] = sort(var(D(:,2:end)));
        D(:,2:end) = D(:,idx+1);
    end
    varOfD(epochs(group), trajbound(group)+1, rewarded(group)+1, leftright(group)+1, :) = var(D(:,2:end));

    if Opt.interpolate ~= 0
        norm = @(x) (x-min(x))./range(x);
        interpfun = @(j) interp1(norm(1:size(D,1)), D(:,j), norm(1:Opt.interpolate),'spline');
        Dnew = arrayfun(interpfun, 2:size(D,2), 'UniformOutput', false);
        Dnew = cellfun(@(x) x', Dnew, 'UniformOutput', false);
        D = [(1:Opt.interpolate)',  cat(2, Dnew{:})];
    end


    %cla
    prev= [];
    for j = 1:size(D,1)
        c = D(j,1);
        d = D(j,dim+1);
        d = num2cell(d);
        s =scatter3(d{:}, 'MarkerFaceColor', colors(c,:), 'MarkerEdgeColor', colors(c,:));
        set(s,'SizeData',10);
        if Opt.distinguishRight && leftright(group)
            set(s,'MarkerFaceAlpha',0.5, 'MarkerEdgeColor', 'k');
        end
        if ~isempty(prev)
            l = num2cell([d{:}; prev{:}],1);
            line(l{:}, 'Color', colors(c,:));
            if Opt.showDirection
                
            end
        end
        hold on;
        prev = d;
    end
    xlim([mx Mx,]);
    ylim([my My,]);
    zlim([mz Mz]);
    title(sprintf('Epoch %d, %d-bnd, \nrew=%d, leftright=%d', epochs(group), trajbound(group), rewarded(group), leftright(group)));
    if Opt.drawnow
        drawnow
    end
end
sgtitle("PSID");
