%%
% This script looks at how animal's performance might be shaping the
% variance of the components' correlation to neural firing

% LIST OF ANIMALS
%%
strength_variance = cell(numel(animal_list),4);
mean_component_strength = cell(numel(animal_list),4);
mu_ci = cell(numel(animal_list),4);

for phase = 1:4
    for anim = 1:numel(animal_list)
        
        beh = table.behavior.lookup(animal_list(anim), animal_critical_behaviors{anim}(phase).component_time);
        comp = animal_critical_behaviors{anim}(phase).component_strength';
        groups = findgroups(beh.epoch);
        
        num_epochs                          = numel(unique(groups));
        strength_variance{anim,phase}       = zeros(1,num_epochs);
        mean_component_strength{anim,phase} = zeros(1,num_epochs);
        
        TIME = 2;
        COMP = 1;
        for u = 1:numel(unique(groups))
            C = comp(:,groups == u);
            T = beh{a,phase}(groups == u, :).tperf_timewise;
            mean_component_strength{anim,phase}(u) = nanmean(nanmean(comp(:, groups == u), TIME), COMP);
            %strength_variance{anim,phase}(u)       = nanstd(mean(comp(:, groups == u),1));
            %strength_variance{anim,phase}(u)       = nanstd(mean(comp(:, groups == u),1));
            strength_variance{anim,phase}(u) = mean(var(bsxfun(@minus, C, median(C,COMP)),[],TIME)); % var( fluctuation of component above its median )
            var_tperf{anim,phase}(u)         = var(bsxfun(@minus, T, median(T,1)),[],1); % var( fluctuation of component above its median )
            mu_tperf{anim,phase}(u)          = mean(T); % var( fluctuation of component above its median )
            abs_diff_tperf{anim,phase}(u)    = mean(abs(diff(T))); % var( fluctuation of component above its median )
            mu_diff_tperf{anim,phase}(u)     = abs(mean(diff(T))); % var( fluctuation of component above its median )
            % RY : I fixed this above : runs

            % RYAN : still need ci for this variance thing above.
            bootfun = @(x) mean(var(bsxfun(@minus, x, median(x,COMP)),[],TIME));
            var_ci{anim,phase}(u,:) = bootci(1000, bootfun, comp(:,groups == u));

        end
    end
end

%% variance and ci of variance
fig('variance of component strength across critical behaviors'); clf; tiledlayout(2,2);
titles = ["reward", "error", "inboundchoice","outboundchoice"];

for phase = 1:4
    nexttile;
    for anim = 1:numel(animal_list)
        plot(strength_variance{anim,phase});
        hold on
    end
    title(titles(phase));
    legend(animal_list)
    xlabel("epoch")
    ylabel("component variance")
end

fig('norm variance of component strength across critical behaviors'); clf; tiledlayout(2,2);
titles = ["reward", "error", "inboundchoice","outboundchoice"];
for phase = 1:4
    nexttile;
    for anim = 1:numel(animal_list)
        S = strength_variance{anim,phase};
        norm = @(x) (x-min(S))./(max(S)-min(S));
        y = norm(S);
        x = 1:numel(y);
        yci = var_ci{anim, phase};
        plot(x,y);
        X = [x;x(end:-1:1)]';
        hold on
    end
    title(titles(phase));
    
    xlabel("epoch")
    ylabel("component variance")
end

%%
fig('norm variance versus norm var tperf'); clf; tiledlayout(2,2);
titles = ["reward", "error", "inboundchoice","outboundchoice"];
for phase = 1:4
    nexttile;
    for anim = 1:numel(animal_list)
        S = strength_variance{anim,phase};
        T = mu_tperf{anim,phase};
        norm = @(x) (x-min(x))./(max(x)-min(x));
        x = norm(T);
        filt = x > 0.20; % throw out early statespace because early tends to be inaccurate, although can see the trend without this
        %filt = true(size(T));
        y = norm(S);
        scatter(x(filt),y(filt),'k')
        %set(fill(X(:), norm(yci(:)),'k'), 'facealpha', 0.2);
        hold on
    end
    title(titles(phase));
    
    xlabel("state space performance" + newline + "(0\rightarrow minimum," + newline+ "1\rightarrow maximum)")
    ylabel("component variance")
    sgtitle("high/stable performance potentially associated with " + newline +  "lower movement in along commm subspace dim.")
end
%

fig('color by animal -norm variance versus norm var tperf'); clf; tiledlayout(2,2);
titles = ["reward", "error", "inboundchoice","outboundchoice"];
colors = crameri('roma', numel(animal_list));
for phase = 1:4
    nexttile;
    for anim = 1:numel(animal_list)
        S = strength_variance{anim,phase};
        T = mu_tperf{anim,phase};
        norm = @(x) (x-min(x))./(max(x)-min(x));
        x = norm(T);
        filt = x > 0.10; % throw out early statespace because early tends to be inaccurate - totally visiable without this, but taking times were statespace has had multiple samples and where the animal has experienced errors to know there's something to learn.
        %filt = true(size(T)); 
        y = norm(S);
        s = scatter(x(filt),y(filt));
        s.MarkerFaceColor = colors(anim, :);
        s.MarkerEdgeColor = colors(anim, :)/2;
        %set(fill(X(:), norm(yci(:)),'k'), 'facealpha', 0.2);
        hold on
    end
    title(titles(phase));
    
    xlabel("state space performance" + newline + "(0\rightarrow minimum," + newline+ "1\rightarrow maximum)")
    ylabel("component variance")
    sgtitle("high/stable performance potentially associated with " + newline +  "lower movement in along commm subspace dim.")
end
%



%% ---
%% CIs
%% ---
fig('Confidence interval across critical behaviors'); clf; tiledlayout(2,2);
HIGH = 2;
LOW = 1;
for phase = 1:4
    nexttile;
    for anim = 1:numel(mu_ci{phase})
        hold on
        plot(var_ci{phase}{anim}(HIGH)-mu_ci{phase}{anim}(LOW));
    end
    
    title(titles(phase));
    xlabel("time")
    ylabel("component strength ci")

end
