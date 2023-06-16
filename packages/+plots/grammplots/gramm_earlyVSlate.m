% this script plots the early and late epochs' comparison
% patterns should first be reshaped into Patterns_AllAnimals
%%
T = query.getPatternTable(Patterns_AllAnimals);

patternType = T.patternType;
patternType( ~contains(patternType, '-' ), : ) = patternType( ~contains(patternType, '-' ), : ) + "-";
columns = patternType.split('-');
columns(columns(:,2)=="",2) = "pattern activity";
T.patternAbstract = columns(:,1);
T.control = columns(:,2);
T.patternAbstractSymbol = T.patternAbstract.replace('theta','$\theta$').replace('delta','$\delta$').replace('ripple','SPW-R');
% T.winSize = cell2mat(T.winSize);


%%  perc dim
perc_dim_skewness = zeros(2,3);
text_size = 8;
point_size = 1.5;

clf
fig('early - late epoch comparison, perc dimension')
clear g

hpcearlysubset = T.epoch == 'early'...
    & T.directionality == "hpc-hpc";
pfcearlysubset = T.epoch == 'early'...
    & T.directionality == "hpc-pfc";
hpclatesubset = T.epoch == 'late'...
    & T.directionality == "hpc-hpc";
pfclatesubset = T.epoch == 'late'...
    & T.directionality == "hpc-pfc";


hpcearlysubset = T(hpcearlysubset,:);
pfcearlysubset = T(pfcearlysubset,:);

hpclatesubset = T(hpclatesubset,:);
pfclatesubset  = T(pfclatesubset,:);


x = hpcearlysubset.percMax_rrDim;
y = hpclatesubset.percMax_rrDim;
perc_dim_skewness(1,1) = skewness(y-x);

corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,1) = gramm( 'x', x,...
    'y', y,...
    'color', categorical(hpcearlysubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(1,1).facet_grid(categorical(hpcearlysubset.control), categorical(hpcearlysubset.patternAbstract))
g(1,1).geom_point('dodge', 0.5, 'alpha', 0.3);
g(1,1).geom_abline('style','k:');
g(1,1).stat_cornerhist(corner_kws{:});
g(1,1).set_point_options('base_size', point_size);
g(1,1).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(1,1).set_names('x', "Early Epochs" + newline + "Pred Dims", ...
    'y', "Late Epochs" + newline + "Pred Dims", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();

y = pfcearlysubset.percMax_rrDim;
x = pfclatesubset.percMax_rrDim;
perc_dim_skewness(1,2) = skewness(y-x);
corner_kws = {'edges',-0.5:0.05:0.5, 'aspect',0.6, 'location', [0.5], 'fill', 'transparent', 'normalization',  'countdensity'}; %WARNING : EDGES OF HISTOGRAM ARE SET MANUALLY! you need this to surround your data
g(1,2) = gramm( ... 'subset', subset,...
    'x', x,...
    'y', y,...
    'color', categorical(pfcearlysubset.patternAbstract));
% assert(all(hpcsubset.patternType == pfcsubset.patternAbstract))
g(1,2).facet_grid(categorical(pfcearlysubset.control), categorical(pfcearlysubset.patternAbstract))
g(1,2).geom_point('dodge', 0.5, 'alpha', 0.3);
g(1,2).geom_abline('style','k:');
g(1,2).stat_cornerhist(corner_kws{:});
g(1,2).set_point_options('base_size', point_size);
g(1,2).set_text_options('label_scaling', 1.2, 'base_size', text_size);

g(1,2).set_names('x', "Early Epochs" + newline + "Pred Dims", ...
    'y', "Late Epochs" + newline + "Pred Dims", ...
    'Color', 'Pattern', ...
    'Lightness', 'Treatment/Control')
g.axe_property('PlotBoxAspectRatio', [1,1,1], 'DataAspectRatioMode','auto')
g.draw();
%% stats
patternPercSum = zeros(2,2,nPatterns);
all_perc = zeros(2,2,3,300);

for m = 1:2 % early vs late
    
    for d = 1:2
        for i = 1:3
            for p = 1:300
                curr = Patterns_AllAnimals(m,p,d,i);
                nTarget_temp = numel(curr.index_target);
                nSource_temp = numel(curr.index_source);
                nDimMax = min(nTarget_temp, nSource_temp);
                currDim = curr.rankRegress.optDimReducedRankRegress;
                currPerc = currDim/nDimMax;
                
                %                 all_theta(m,d,:) = [all_theta(m,d,:), currPerc];
                patternPercSum(m,d,i) = patternPercSum(m,d,i) ...
                    + currPerc;
                all_perc(m,d,i,p) = currPerc;
                
            end
        end
    end
end

patternPercSum = patternPercSum./300;
all_perc_mean = mean(all_perc, 4);
all_perc_std = std(all_perc,0,4);

for i = 1:3
    for j = 1:2
        [h_epoch(i,j), p_epoch(i,j)] = ttest2 (all_perc(1,j,i,:),all_perc(2,j,i,:));
        mean_early(j,i) = mean(all_perc(1,j,i,:));
        mean_late(j,i) = mean(all_perc(2,j,i,:));
        
        std_early(j,i) = std(all_perc(1,j,i,:));
        std_late(j,i) = std(all_perc(1,j,i,:));
    end
end
%%
fig("comparing early vs. late dimensions"), clf

% bar plots of all rhythm mean across animals and std
x = 1:3;

bar(x,[mean_early(1,:); mean_late(1,:)])
alpha(0.33)
hold on
er = errorbar(x,mean_early(1,:), std_early(1,:));
er2 = errorbar(x,mean_late(1,:), std_late(1,:));

xticklabels(patternnames)
ylabel("sum of perc. dims spanned")

legend("early","late")


legend("early","late")
%%
fig("comparing early vs. late dimensions pfc"), clf

% bar plots of all rhythm mean across animals and std
x = 1:3;

bar(x,[mean_early(2,:); mean_late(2,:)])
alpha(0.33)
hold on
er = errorbar(x,mean_early(2,:), std_early(2,:));
er2 = errorbar(x,mean_late(2,:), std_late(2,:));

xticklabels(patternnames)
ylabel("sum of perc. dims spanned")

legend("early","late")

