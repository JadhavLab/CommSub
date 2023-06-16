%%
fig('variance of all comp strength across epochs '+animal);
clf; tiledlayout(2,2);
titles_bypattern = ["all","theta","delta","ripple"];
for i = 1:4
    
    nexttile
    
    for j = 1:numel(mean_hpc_variance{i})
        hold on
        stem(j:j,mean_hpc_variance{i}(j),"filled","color", "black");
        hold on
        stem(j:j,mean_pfc_variance{i}(j),"color", "black");
        
    end
    title(titles_bypattern(i))
    hold on
    legend("hpc-hpc","hpc-pfc")
    ylabel("variance")
    xlabel("epoch")
    
end
%%
fig('variance of comp strength across epochs, by dim '+animal);
clf; tiledlayout(2,2);
titles_bydim = ["dim 1","dim 2","dim 3"];
for i = 8:10
    
    nexttile
    
    for j = 1:numel(mean_hpc_variance{i})
        hold on
        stem(j:j,mean_hpc_variance{i}(j),"filled","color", "black");
        hold on
        stem(j:j,mean_pfc_variance{i}(j),"color", "black");
        
    end
    title(titles_bydim(i-7))
    hold on
    legend("hpc-hpc","hpc-pfc")
    ylabel("variance")
    xlabel("epoch")
    
end
%%
slope_hpc_sv_mean = mean(slope_hpc_sv_all, 2);
slope_pfc_sv_mean = mean(slope_pfc_sv_all, 2);
slope_hpc_sv_std = std (slope_hpc_sv_all,0, 2);
slope_pfc_sv_std  = std (slope_pfc_sv_all,0, 2);
%%
fig('variance change fit across epochs, by pattern, hpc '+animal); clf

bar( slope_hpc_sv_mean(1:4));


fig('variance change fit across epochs, by pattern, pfc '+animal); clf

errorbar(slope_pfc_sv_mean(1:4),slope_pfc_sv_std(1:4));

%% old-partitioned analysis on phases of behaviors

%% variance by phase
for p = 1: Option.numPartition
    for phase = 1:4
        
        
        hpcbeh = table.behavior.lookup(animal, Patterns(p,1,1).compAnalysisByBeh(phase).component_time);
        
        pfcbeh = table.behavior.lookup(animal, Patterns(p,2,1).compAnalysisByBeh(phase).component_time);
        
        hpccomp = Patterns(p,1,1).compAnalysisByBeh(phase).component_strength';
        pfccomp = Patterns(p,2,1).compAnalysisByBeh(phase).component_strength';
        groups = findgroups(hpcbeh.epoch);
        
        num_epochs                          = numel(unique(groups));
        hpc_strength_variance{p}{phase}       = zeros(1,num_epochs);
        pfc_strength_variance{p}{phase}       = zeros(1,num_epochs);
        % hpc_mean_component_strength{phase} = zeros(1,num_epochs);
        
        TIME = 2;
        COMP = 1;
        for u = 1:numel(unique(groups))
            C_hpc = hpccomp(:,groups == u);
            C_pfc = pfccomp(:,groups == u);
            T_hpc = hpcbeh(groups == u, :).tperf_timewise;
            T_pfc = pfcbeh(groups == u, :).tperf_timewise;
            pfc_strength_variance{p}{phase}(u)       = mean(var(bsxfun(@minus, C_pfc, median(C_pfc,COMP)),[],TIME));
            
            hpc_strength_variance{p}{phase}(u)       = mean(var(bsxfun(@minus, C_hpc, median(C_hpc,COMP)),[],TIME));
            
        end
    end
    
end

%%
for p = 1:Option.numPartition
    for i = 1:4
        hpc_variance{i}(p,:) = hpc_strength_variance{p}{i};
        pfc_variance{i}(p,:) = pfc_strength_variance{p}{i};
    end
end

for i = 1:4
    mean_hpc_variance{i} = mean(hpc_variance{i},1);
    mean_pfc_variance{i} = mean(pfc_variance{i},1);
    
    std_hpc_variance{i} = std(hpc_variance{i});
    std_pfc_variance{i} = std(pfc_variance{i});
end

%%
fig('variance of comp strength across epochs '+animal);
clf; tiledlayout(2,2);

for i = 1:4
    
    nexttile
    
    for j = 1:numel(mean_hpc_variance{i})
        hold on
        stem(j:j,mean_hpc_variance{i}(j),"filled","color", "black");
        hold on
        stem(j:j,mean_pfc_variance{i}(j),"color", "black");
        
    end
    
    hold on
    legend("hpc-hpc","hpc-pfc")
    ylabel("variance")
    xlabel("epoch")
    
end