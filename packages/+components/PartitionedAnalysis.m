%% across all times, component strength
for p = 1: Option.numPartition
    hpcbeh = table.behavior.lookup(animal, Patterns(p,1,1).compAnalysisTable.time);
    
    hpccomp = Patterns(p,1,1).compAnalysisTable.subspaceComponents';
    pfccomp = Patterns(p,2,1).compAnalysisTable.subspaceComponents';
    groups = findgroups(hpcbeh.epoch);
    
    num_epochs                          = numel(unique(groups));
    for k = 1:12 % 1: all, 2:theta, 3:delta, 4:ripple
        hpc_strength_variance_all{p}{k}      = zeros(1,num_epochs);
        pfc_strength_variance_all{p}{k}      = zeros(1,num_epochs);
    end
    % hpc_mean_component_strength{phase} = zeros(1,num_epochs);
    
    TIME = 2;
    COMP = 1;
    for u = 1:numel(unique(groups))
        C_hpc = hpccomp(:,groups == u);
        C_pfc = pfccomp(:,groups == u);
        %
        %         T_hpc = hpcbeh(groups == u, :).tperf_timewise;
        %         T_pfc = pfcbeh(groups == u, :).tperf_timewise;
        % mean_component_strength{anim,phase}(u) = nanmean(nanmean(comp(:, groups == u), TIME), COMP);
        %strength_variance{anim,phase}(u)       = nanstd(mean(comp(:, groups == u),1));
        %strength_variance{anim,phase}(u)       = nanstd(mean(comp(:, groups == u),1));
        pfc_strength_variance_all{p}{1}(u)       = mean(var(bsxfun(@minus, C_pfc, median(C_pfc,COMP)),[],TIME)); % var( fluctuation of component above its median )
        
        hpc_strength_variance_all{p}{1}(u)       = mean(var(bsxfun(@minus, C_hpc, median(C_hpc,COMP)),[],TIME)); % var( fluctuation of component above its median )
        
        for k = 2:7
            C_hpc_pat = C_hpc((k-2)*5+1:(k-1)*5,:);
            C_pfc_pat = C_pfc((k-2)*5+1:(k-1)*5,:);
            pfc_strength_variance_all{p}{k}(u)       = mean(var(bsxfun(@minus, C_pfc_pat, median(C_pfc_pat,COMP)),[],TIME)); % var( fluctuation of component above its median )
            hpc_strength_variance_all{p}{k}(u)       = mean(var(bsxfun(@minus, C_hpc_pat, median(C_hpc_pat,COMP)),[],TIME)); % var( fluctuation of component above its median )
            pfc_m_all(p,k,u) = nanmean(nanmean(C_pfc_pat, TIME), COMP);
            hpc_m_all(p,k,u) = nanmean(nanmean(C_hpc_pat, TIME), COMP);
        end
        
        
        
        for k = 8:12
            C_hpc_dim = C_hpc((k-8)*6+1:(k-7)*6,:);
            C_pfc_dim = C_pfc((k-8)*6+1:(k-7)*6,:);
            pfc_strength_variance_all{p}{k}(u)    = mean(var(bsxfun(@minus, C_pfc_dim, median(C_pfc_dim,COMP)),[],TIME));
            hpc_strength_variance_all{p}{k}(u) = mean(var(bsxfun(@minus, C_hpc_dim, median(C_hpc_dim,COMP)),[],TIME));
            pfc_m_all(p,k,u) = nanmean(nanmean(C_pfc_dim, TIME), COMP);
            hpc_m_all(p,k,u) = nanmean(nanmean(C_hpc_dim, TIME), COMP);
        end
    end
    
end
% make table for mean strength by epoch


%% organize results
for p = 1:Option.numPartition
    for i = 1:12
        hpc_variance{i}(p,:) = hpc_strength_variance_all{p}{i};
        pfc_variance{i}(p,:) = pfc_strength_variance_all{p}{i};
    end
end

for i = 1:12
    mean_hpc_variance{i} = mean(hpc_variance{i},1);
    mean_pfc_variance{i} = mean(pfc_variance{i},1);
    
    std_hpc_variance{i} = std(hpc_variance{i});
    std_pfc_variance{i} = std(pfc_variance{i});
end


%% poly-fits for changes in component strength variance

for j = 1:12
    for p = 1:Option.numPartition
        P_hpc_sv_all{j}{p} = polyfit(1:numel(hpc_variance{j}(p,:)),hpc_variance{j}(p,:),1);
        P_pfc_sv_all{j}{p} = polyfit(1:numel(pfc_variance{i}(p,:)),pfc_variance{j}(p,:),1);
        
        slope_hpc_sv_all(j,p) = P_hpc_sv_all{j}{p}(1);
        slope_pfc_sv_all(j,p) = P_pfc_sv_all{j}{p}(1);
    end
end

%store into a table
animalName = repmat(animal, 2*Option.numPartition, 1);
inxName = [repmat("hpc-hpc", Option.numPartition, 1);repmat("hpc-pfc", Option.numPartition, 1)];
SlopeTable = array2table([animalName, inxName, [slope_hpc_sv_all' ; slope_pfc_sv_all']]);
SlopeTable.Properties.VariableNames(1:14) = ...
    {'animal','inx','all','theta','delta','ripple','low-theta','low-delta','low-ripple','dim1','dim2','dim3','dim4','dim5'};
load('PolyFitTable_AllTimes.mat');
PolyFitTable_AllTimes = [PolyFitTable_AllTimes; SlopeTable];
save('PolyFitTable_AllTimes','PolyFitTable_AllTimes','-v7.3')
% Default heading for the columns will be A1, A2 and so on.
% You can assign the specific headings to your table in the following manner

%% component strength mean strength during different rhythms

hpc_m_allepochs = mean(hpc_m_all,3);
pfc_m_allepochs = mean(pfc_m_all,3);

hpc_meanCpStrength = mean(hpc_m_allepochs,1);
pfc_meanCpStrength = mean(pfc_m_allepochs,1);

hpc_stdCpStrength = std(hpc_m_allepochs,0,1);
pfc_stdCpStrength = std(pfc_m_allepochs,0,1);


%%
% reorganize the matrix in different ways and then plot imagesc of the
% correlation

for a = 1:size(Behaviors,1)
    
    
    for p = 1:50
        for d = 1:2
            all_subspaces = Behaviors(a,p,d).allSubSpaceComponents;
            subspace_by_dim = components.organizeComponents(all_subspaces, 5, 3, "dim");
            for i = 1:3
                dim_pattern = subspace_by_dim((i-1)*6+1:(i-1)*6+3,:);
                dim_corr = abs(corrcoef(dim_pattern'));
                all_corr(a,p,d,i) = sum(dim_corr,'all');
            end
            
            for m = 1:4
                
                components_by_dim = components.organizeComponents...
                    (Behaviors(a,p,d).compAnalysisByBeh(m).critical_components', 5, 3, "dim");
                
                for k = 1:3
                    dim_pattern_bybeh = components_by_dim((k-1)*6+1:(k-1)*6+3,:);
                    dim_pattern_bybeh(isnan(dim_pattern_bybeh)) = 0;
                    dim_corr_bybeh = abs(corrcoef(dim_pattern_bybeh'));
                    
                    all_corr_byBeh(a,p,d,k,m) = sum(dim_corr_bybeh,'all');
                end

            end
        end
    end
end

%%
%%
all_corr_sumed = squeeze(mean(squeeze(mean(all_corr,2)),1));

all_corr_sumed_bybeh =  squeeze(mean(squeeze(mean(all_corr_byBeh,2)),1));


allpat_corr = reshape(all_corr, [300,2,3]);
allpat_corr_bybeh = reshape(all_corr_byBeh, [300,2,3,4]);

allpat_corr_bybeh(:,:,:,5) = allpat_corr;

hpccorr = allpat_corr_bybeh(:,1,:,:);
pfccorr = allpat_corr_bybeh(:,2,:,:);
hpccorr = squeeze(hpccorr);
pfccorr = squeeze(pfccorr);



for i = 1:3
    for j = 1:5
    mean_corr_hpc(i,j) = mean(hpccorr(:,i,j));
    std_corr_hpc(i,j) = std(hpccorr(:,i,j));
    
        mean_corr_pfc(i,j) = mean(pfccorr(:,i,j));
    std_corr_pfc(i,j) = std(pfccorr(:,i,j));
    end
end

fig('correlation by beh and for all times')

for i = 1:3
    subplot(1,3,i)
    bar(1:5, mean_corr_hpc(i,:))
    hold on
    bar(1:5, mean_corr_pfc(i,:))
end
%%
for i = 1:4
    for p = 1:Option.numPartition
        temp = R_bybeh{p}{i}(1:6,1:6);
        avg_leadingDimCorr(p,:) = mean(abs(temp), 1);
    end
    beh_leadingDimCorr(i,:) = mean(avg_leadingDimCorr,1);
end


subspace_by_dim = components.organizeComponents(behavior_running_subspaces, 5, 3, "dim");
R = abs(corrcoef(subspace_by_dim'));

mean(R(1:6,1:6),1)
