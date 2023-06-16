function plots(psid_struct, beh, varargin)
% Runs a battery of basic plots on psid results, per epoch

ip = inputParser;
ip.addParameter('beh_fields', []);
ip.addParameter('traj_smoothing', true);
ip.addParameter('overall_projLowerD', 3);
ip.addParameter('correctOnly', false);
ip.parse(varargin{:})
Opt = ip.Results;

beh_fields = Opt.beh_fields;
epochs = beh.epoch;
traj = beh.traj;
epochGroups =  findgroups(epochs);
uepochGroups = unique(epochGroups);

for field = string(fieldnames(psid_struct))'
    assign(field, psid_struct.(field));
end
if ~isfield(psid_struct, "set_name")
    set_name = ""
else
    set_name = set_name + " ";
end
clear psid_struct;

keyboard;
if Opt.correctOnly
    behNew = [];
    for g = uepochGroups'
        gBeh = beh(epochGroups == g,:);
        correct =  gBeh.rewarded == 1;
        Z_test{g}(~correct) = [];
        Z_pred{g}(~correct) = [];
        X_pred{g}(~correct) = [];
        gBeh(~correct,:) = [];
        behNew = [behNew; gBeh];
    end
    beh = behNew;
end


%
%fig(set_name + 'Correlation of actual to pred');
%clf
%tiledlayout('flow')
%for g = uepochGroups'
%    nexttile
%    zz = [Z_test{g},  Z_pred{g}];
%    cZ = corrcoef(zz);
%    labels = [["actual_{" + beh_fields + "}"], ["pred_{" + beh_fields + "}"]];
%    heatmap(labels,labels,abs(cZ));
%    crameri('tokyo','pivotValue',0.5);
%    title("Epoch " + g*2);
%end
%
%%% Plot results : Traj
%%% -------------------
%fig(set_name + 'Trajectory');
%ex = 3;
%plot_dims = {ex, numel(uepochGroups)};
%tile = tiledlayout(plot_dims{:})
%set(tile, 'TileIndexing', 'ColumnMajor', 'TileSpacing', 'Tight');
%for g = uepochGroups'
%    trajgroups = findgroups(traj(epochGroups == g));
%    bg = beh(epochGroups == g, :);
%    utrajgroups = unique(trajgroups);
%    epoch_z =  X_pred{g};
%    tcnt = 0;
%    for t = utrajgroups(randperm(numel(utrajgroups)))'
%        tcnt = tcnt+1;
%        if tcnt > ex
%            continue
%        end
%        epoch_traj_z = epoch_z(trajgroups == t,:);
%        epoch_traj_z = num2cell(epoch_traj_z, 1);
%        nexttile;
%        if Opt.traj_smoothing
%            epoch_traj_z = cellfun(@smoothdata, epoch_traj_z, 'UniformOutput', false);
%        end
%        plot3(epoch_traj_z{end-2:end}, ':o', 'MarkerSize', 4);
%        title("Rewarded = " + bg(trajgroups == t,:).rewarded(1) + newline + "Trajbound = " + bg(trajgroups == t,:).trajbound(1) ...
%            + newline + "Left = " + bg(trajgroups==t,:).leftright(1))
%    end
%end
%minmax = @(x) [min(x(:,1)), max(x(:,2))];
%axs = findobj(gcf,'type','axes');
%xlim = get(axs, 'xlim');
%xlim = cat(1, xlim{:});
%ylim = get(axs, 'ylim');
%ylim = cat(1, ylim{:});
%zlim = get(axs, 'zlim');
%zlim = cat(1, zlim{:});
%set(axs, 'xlim', minmax(xlim), 'ylim', minmax(ylim), 'zlim', minmax(zlim));
%
%
%%% Plot results : Overall
%%% ----------------------
%for field = beh_fields
%    fig(set_name + 'Overall' + " Z: Behavioral Projection" + newline + field)
%    clf
%    tiledlayout('flow')
%    norm = @(x) (x-min(x))./(max(x)-min(x));
%    for g = uepochGroups'
%        nexttile
%        zz = bsxfun(@rdivide, bsxfun(@minus, Z_pred{g}, mean(Z_pred{g})), range(Z_pred{g}));
%        [~,S,V] = svd(zz' * zz);
%        if Opt.overall_projLowerD
%            zz = projOrth(zz, V(1:Opt.overall_projLowerD,:));
%        end
%        zz = num2cell(zz,1);
%        labels = beh_fields;
%        lset = 1:3;
%        S = 15*ones(size(beh.(field)(epochGroups==g)));
%        C = norm(beh.(field)(epochGroups==g));cla
%        pivot = 0.5;
%        try
%            colormap(psid.behaviorColor(field, 'pivotValue', pivot));
%        catch
%            colormap(crameri('bukavu', 'pivotValue', pivot));
%        end
%        scatter3(zz{lset}, S, C, 'filled');
%        set(gca, 'Color', 0.3 * ones(1,3))
%        alpha(0.6);
%        xlabel(lset(1));
%        ylabel(lset(2));
%        zlabel(lset(3));
%        title("Epoch " + g*2);
%    end
%end


%% Plot results : Overall
%% ----------------------
for field = beh_fields
    fig(set_name + 'Field Overall'+ " X: Internal State" + newline + field)
    clf
    tiledlayout('flow')
    norm = @(x) (x-min(x))./(max(x)-min(x));
    for g = uepochGroups'
        nexttile
        keyboard
        xx = bsxfun(@rdivide, bsxfun(@minus, X_pred{g}, mean(X_pred{g})), range(X_pred{g}));
        [~,S,V] = svd(xx' * xx);
        if Opt.overall_projLowerD
            xx = projOrth(xx, V(1:Opt.overall_projLowerD,:));
        end
        xx = num2cell(xx,1);
        labels = beh_fields;
        lset = 1:3;
        try
            colormap(psid.behaviorColor(field, 'pivotValue', 0.5));
        catch
            colormap(crameri('bukavu', 'pivotValue', 0.5));
        end
        S = 15*ones(size(beh.directional_lindist_imag(epochGroups==g)));
        C = norm(beh.(field)(epochGroups==g));cla
        scatter3(xx{lset}, S, C, 'filled');
        set(gca, 'Color', 0.3 * ones(1,3))
    %     keyboard
        alpha(0.6);
        xlabel(lset(1));
        ylabel(lset(2));
        zlabel(lset(3));
        title("Epoch " + g*2);
    end
end

%% Plot results : Overall
%% ----------------------
for g = uepochGroups'
    fig(set_name + 'Epoch Overall'+ " X: Internal State" + newline + "Epoch : " + g*2)
    clf
    tiledlayout('flow')
    norm = @(x) (x-min(x))./(max(x)-min(x));
    for field = beh_fields
        nexttile
        xx = bsxfun(@rdivide, bsxfun(@minus, X_pred{g}, mean(X_pred{g})), range(X_pred{g}));
        [~,S,V] = svd(xx' * xx);
        if Opt.overall_projLowerD
            xx = projOrth(xx, V(1:Opt.overall_projLowerD,:));
        end
        xx = num2cell(xx,1);
        labels = beh_fields;
        lset = 1:3;
        try
            colormap(psid.behaviorColor(field, 'pivotValue', 0.5));
        catch
            colormap(crameri('bukavu', 'pivotValue', 0.5));
        end
        S = 15*ones(size(beh.directional_lindist_imag(epochGroups==g)));
        C = norm(beh.(field)(epochGroups==g));cla
        scatter3(xx{lset}, S, C, 'filled');
        set(gca, 'Color', 0.3 * ones(1,3))
    %     keyboard
        alpha(0.6);
        xlabel(lset(1));
        ylabel(lset(2));
        zlabel(lset(3));
        title("Behavior: " + field);
    end
end
