% paths
addpath(genpath(codedefine()))
addpath(hashdefine())

%% load
multi_epoch = false; % usually first 3, last 3 epochs
load("RunsSummary.mat");
% TABLE : ACQUIRE RUNS
% ----------------------
% Determine keys to use : you can use this string to arbitrarily select rows
%   each item of the filtstring is a property to select. $x pulls the x
%   column and applies the test shown
filtstring = ["ismember($animal, [""JS21"",""ZT2"",""ER1"",""JS14"",""JS13"",""JS17""])",...
                                   ..."$spikeBinSize==0.15",...
                                   "$numPartition==50",...
                                   "$quantileToMakeWindows == 0.85",...
                                   "arrayfun(@(x)isequal($winSize(x,:), [0,0.3]), 1:size($winSize,1))'"];
% Get the proper keys
matching_runs = query.getHashed_stringFilt(RunsSummary, filtstring);
disp("Number of matches found: " + height(matching_runs))

%%

%%

if multi_epoch
    keys = [];
    for i = 1:height(matching_runs)
        if numel(matching_runs(i,:).epochsSelected{1} ) ==3 % change this to 3 for 3 early and 3 late
            keys = [keys, matching_runs(i,:).hash];
        end
    end
else
    keys = matching_runs.hash;
end

%% Load keyed save files from local or server
% Load up combined pattern data
localmode = true; % set to false to load from server
if localmode
    % Load locally
    Out = query.combinePatterns(keys, ...
    'RunsSummary', RunsSummary);
    Patterns = Out.Patterns;
    Option = Out.Option;
    clear Out
else
    % Grab from server and load
    [Patterns, otherData] = query.pattern.tmpLoadAndCombine(keys);
    Option = otherData{1}.Option;
end
assert(~isempty(Option), "Data is empty -- downstream code will fail")
Patterns = nd.merge(Patterns, Option, 'only', {'animal', 'generateH', 'genH_name'}, ...
                                                'broadcastLike', true,...
                                                 'overwrite', true);
[nDataset, nPartition, ~, nResult] = size(Patterns)
nPatterns = nResult/2;

uAnimals = unique([Option.animal]);

nSource = zeros(1,nDataset);
nTarget = zeros(1,nDataset);
numDimsUsedForPrediction = cell(1,nDataset);
for a = 1:nDataset
    nSource(a) = size(Patterns(a,1,1,1,1).X_source,1);
    nTarget(a) = size(Patterns(a,1,1,1,1).X_target,1);
    numDimsUsedForPrediction{a} = 1:min(nSource(a), nTarget(a));
end

Patterns = util.type.castefficient(Patterns, ...
            'compressReals', true, 'verbose', false);
Patterns = nd.apply(Patterns, "nSource-X_source", @(x) size(x,1) );
Patterns = nd.apply(Patterns, "nTarget-X_target", @(x) size(x,1) );

% little check
[[Patterns(:, 1,1,1,1,1,1).animal]', [Option(:,1,1,1,1,1).animal]']
[[Patterns(:, 1,1,1,1,1,1).generateH]', [Option(:,1,1,1,1,1).generateH]']
P = munge.reshapeByLabels(Patterns, 1, [Option.generateH],  'checksumSplitField', 'animal');
O   = munge.reshapeByLabels(Option, 1,   [Option.generateH], 'checksumSplitField', 'animal');


% Patterns = nd.flexrmfield(Patterns, {'X_source', 'X_target'});
dump = matfile(fullfile(codedefine,"figures","SPF"), "Writable", true);

%%  
% The different animals loaded will actually collapse into the partition
% 
% dimenion, segregated by the genH methods
% 
% *(nMethods * nPartiton * nDirection * nPatterns)*

%% 
% ISSUE: what is `Var16`?
T = query.getPatternTable(Patterns);

% Figure 2A: Cofiring in source and target
% 
% Figure 2B: Prediction Performance
% 
% Figure 4A/B: Results from rank regress
% 
% Figure 4C: How optimal number of predicitive dimensions compare between hpc 
% and pfc
% 
% Figure 6: removing predicitive dimensions from both sources
% 
% Figure 7: Predictive performance from dominant dimensions
%% Figure 2
% 
% A: How pairs of neuron in source/target co-fire
% This is method insensitive. Just lump over animals

%% 2A calculate co-firing across all animals
Fig2 = struct();
Fig2.a.all = plots.cf.cofiring(Patterns, Option, 'appendAttributes', {'animal', 'generateH'});
[spec, opt] = munge.getH(Patterns, Option, "spec");
Fig2.a.spec = plots.cf.cofiring(spec, opt, 'appendAttributes', {'animal', 'generateH'});
[coh, opt] = munge.getH(Patterns, Option, "coh");
Fig2.a.coh = plots.cf.cofiring(coh, opt);
[wpli, opt] = munge.getH(Patterns, Option, "wpli");
Fig2.a.wpli = plots.cf.cofiring(wpli, opt, 'appendAttributes', {'animal', 'generateH'});
dump.Fig2 = Fig2;

% ----------------
% Fields created:
% ----------------
%.withhpc_pairs            ; % literally list of correlations between source and target=hpc
%.withpfc_pairs            ; % literally list of correlations between source and target=pfc, but separable for each batch
%.all_pairs_withhpc        ; % same as above, but all in one vector
%.all_pairs_withpfc        ; % same as above, but all in one vector
%.mean_corrwithhpc         ; % mean of all_pairs_withhpc
%.mean_corrwithpfc         ; % mean of all_pairs_withpfc
%.std_corrwithhpc          ; % std of all_pairs_withhpc
%.std_corrwithpfc          ; % std of all_pairs_withpfc
%.mean_withhpccorr_pattern ; % mean of withhpc_pairs
%.mean_withpfccorr_pattern ; % mean of withpfc_pairs
%.prophpc                  ; % properties of hpc

%% 
% The co-firing of hpc and pfc neurons during different activity patterns, with 
% average plotted

% TODO:
% 1. Add sig to figs
% 2. Redo with different normalization
% 3. Check standard coherence w/ fresh chronux
% 4. Show ranges with quantile plot

Fig2.a.all  = plots.cf.plotCofiring(Fig2.a.all, Option, 'figAppend', 'all');
Fig2.a.spec = plots.cf.plotCofiring(Fig2.a.spec,Option, 'figAppend', 'spec');
Fig2.a.coh  = plots.cf.plotCofiring(Fig2.a.coh, Option, 'figAppend', 'coh');
Fig2.a.wp   = plots.cf.plotCofiring(Fig2.a.wpli,Option, 'figAppend', 'wp');

plots.cf.plotCofiring(Fig2.a.all, Option, 'Normalization', 'probability', 'figAppend', 'all');
plots.cf.plotCofiring(Fig2.a.spec,Option, 'Normalization', 'probability', 'figAppend', 'spec');
plots.cf.plotCofiring(Fig2.a.coh, Option, 'Normalization', 'probability', 'figAppend', 'coh');
plots.cf.plotCofiring(Fig2.a.wpli,Option, 'Normalization', 'probability', 'figAppend', 'wp');

plots.cf.plotCofiring(Fig2.a.all, Option, 'Normalization', 'cdf', 'figAppend', 'all');
plots.cf.plotCofiring(Fig2.a.spec,Option, 'Normalization', 'cdf', 'figAppend', 'spec');
plots.cf.plotCofiring(Fig2.a.coh, Option, 'Normalization', 'cdf', 'figAppend', 'coh');
plots.cf.plotCofiring(Fig2.a.wpli,Option, 'Normalization', 'cdf', 'figAppend', 'wp');


%% 
% Example 
% the difference in cofiring between hpc-hpc pairs and hpc-pfc pairs
% are significant for all the activity patterns

plots.cf.plotCfExampByDirection(Fig2.a.all,  Patterns, Option, "figAppend", 'all');
plots.cf.plotCfExampByDirection(Fig2.a.spec, Patterns, Option, "figAppend", 'spec');
plots.cf.plotCfExampByDirection(Fig2.a.coh,  Patterns, Option, "figAppend", 'coh');
plots.cf.plotCfExampByDirection(Fig2.a.wpli, Patterns, Option, "figAppend", 'wp');

plots.cf.plotCfExampByDirection(Fig2.a.all,  Patterns, Option, "figAppend", 'all',  'Normalization', 'probability');
plots.cf.plotCfExampByDirection(Fig2.a.spec, Patterns, Option, "figAppend", 'spec', 'Normalization', 'probability');
plots.cf.plotCfExampByDirection(Fig2.a.coh,  Patterns, Option, "figAppend", 'coh',  'Normalization', 'probability');
plots.cf.plotCfExampByDirection(Fig2.a.wpli, Patterns, Option, "figAppend", 'wp',   'Normalization', 'probability');

%%

dump.Fig2 = Fig2;

%% 
% "These weak correlations indicate that only a small fraction of a neuron’s 
% response variability can be explained by another individual neuron"
% 
% B: Explained Variance
Fig2.b = plots.pred.var.regionalexplained(Patterns, Option, ...
                        'appendAttributes', {'animal', 'generateH','name'});
Fig2.b = plots.pred.var.plotexplained(Fig2.b, Option, ...
                        'figAppend', 'all');
Fig2.b = plots.pred.var.plotexplained(Fig2.b, Option, ...
                        'figAppend', 'all-log',...
                        'yscale', 'log');

% Fieldss created:
% ----------------
% .r_withhpc_patterns; % firing prediction with hpc
% .r_withpfc_patterns; % firing prediction with pfc
% .single_pred_with_hpc; % single neuron firing prediction with hpc
% .single_pred_with_pfc; % single neuron firing prediction with pfc
% .patternVarExplained_hpc; % variance explained with hpc
% .patternVarExplained_pfc; % variance explained with pfc

%% 
% Predicition Performance

%%
% Check how dist of corr differs from dist of pred score 
[h_hpc, p_hpc] = kstest2(Fig2.b.meanpred_hpc, Fig2.a.all.mean_withhpccorr_pattern);
[h_pfc, p_pfc] = kstest2(Fig2.b.meanpred_pfc, Fig2.a.all.mean_withpfccorr_pattern);

%%
mean_withhpc = mean(r_square_withhpc(intersect(~isinf(r_square_withhpc), ~isnan(r_square_withhpc))));
std_withhpc  = std(r_square_withhpc(intersect(~isinf(r_square_withhpc), ~isnan(r_square_withhpc))));

mean_withpfc = mean(r_square_withpfc(intersect(~isinf(r_square_withpfc), ~isnan(r_square_withpfc))));
std_withpfc  = std(r_square_withpfc(intersect(~isinf(r_square_withpfc), ~isnan(r_square_withpfc))));

subplot(2,1,1)
ax1 = nexttile;
hist_withpfc = histogram(r_square_withhpc,25);
ylabel("Data Sets")
title ("source predicting HPC targets")
lineObject=line([mean_withhpc,mean_withhpc],[0 max(hist_withpfc.Values)]);
lineObject.LineStyle = ':'; % Make line dotted
lineObject.LineWidth = 2;  % Thicken the line
lineObject.Color = 'black'; % Color it black

lineObject2 = line([median_singlehh,median_singlehh],[0 max(hist_withpfc.Values)]);
lineObject2.LineWidth = 3;
lineObject2.Color = 'blue';
%%

subplot(2,1,2)
ax2 = nexttile;
hist_hp = histogram(r_square_withpfc,25);
% prediction performance from hpc to pfc cells
ylabel("Data Sets")
title ("source predicting PFC targets")
xlabel("Performance")
linkaxes([ax1,ax2],'x');
lineObject=line([mean_withpfc,mean_withpfc],[0 max(hist_hp.Values)]);
lineObject.LineStyle = ':'; % Make line dotted
lineObject.LineWidth = 2;  % Thicken the line
lineObject.Color = 'black'; % Color it black

lineObject2 = line([median_singlehp,median_singlehp],[0 max(hist_hp.Values)]);
lineObject2.LineWidth = 3;  % Thicken the line
lineObject2.Color = 'blue'; % Color it black

%%
figure(601)
clf
for i = 1:nPatterns
    subplot(3,1,i)
    plot(patternPerformance_pfc(i,:));
    hold on
    plot(patternPerformance_hpc(i,:));
    
    legend("pfc","hpc")
    % how to interpret nans
end
%%

% print stats
formatSpec1 = "%s: %0.3f±%0.3f";
disp("Prediction Performance on Average")
sprintf(formatSpec1,Patterns(1,1).directionality,mean_withhpc,std_withhpc)
sprintf(formatSpec1,Patterns(2,1).directionality,mean_withpfc,std_withpfc)
if useSinglePrediction
    disp('single source prediction median')
    formatSpec2 = "%s: %0.5e";
    sprintf(formatSpec2,Patterns(1,1).directionality,median_singlehh)
    sprintf(formatSpec2,Patterns(1,2).directionality,median_singlehp)
end
%% Figure 4

withPredDims;
%%
% histogram of all predicition perf
hpc_theta = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-hpc" & T.patternType == "theta";
hpc_delta = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-hpc" & T.patternType == "delta";
hpc_ripple = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-hpc" & T.patternType == "ripple";


pfc_theta = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-pfc" & T.patternType == "theta";
pfc_delta = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-pfc" & T.patternType == "delta";
pfc_ripple = T.generateH == 'fromFilteredEEG  fromRipTimes'...
    & T.directionality == "hpc-pfc" & T.patternType == "ripple";

hpc_pred = {T(hpc_theta,:).full_model_performance,T(hpc_delta,:).full_model_performance,T(hpc_ripple,:).full_model_performance};
pfc_pred = {T(pfc_theta,:).full_model_performance,T(pfc_delta,:).full_model_performance,T(pfc_ripple,:).full_model_performance};



h_corrdiff = zeros(1,nPatterns);
p_corrdiff = zeros(1,nPatterns);
for i = 1:3
    subplot(3,1,i);
    hold off;
    predperf_hh = histogram(hpc_pred{i});
    set(predperf_hh, 'EdgeColor', 'none', 'FaceAlpha', 0.33);
    hold on
    predperf_hp = histogram(pfc_pred{i});
    set(predperf_hp, 'EdgeColor', 'none', 'FaceAlpha', 0.33);
    
    title(patternnames(i))
    ylabel("data sets")
    
    hold on
    avg_hh=line([mean(hpc_pred{i}),mean(hpc_pred{i})],[0 max(predperf_hh.Values)]);
    avg_hh.LineStyle = ':'; % Make line dotted
    avg_hh.LineWidth = 2;  % Thicken the line
    avg_hh.Color = 'blue';
    avg.hh.DisplayName = "hpc-hpc corfiring mean";
    %
    hold on
    avg_hp=line([mean(pfc_pred{i}),mean(pfc_pred{i})],[0 max(predperf_hp.Values)]);
    avg_hp.LineStyle = ':'; % Make line dotted
    avg_hp.LineWidth = 2;  % Thicken the line
    avg_hp.Color = 'red';
    avg.hh.DisplayName = "hpc-pfc corfiring mean";
    
    if source == "hpc"
        legend("hpc-hpc","hpc-pfc")
    else
        legend("pfc-hpc","pfc-pfc")
    end
    
    xlabel("predicition performance")
    [h_corrdiff(i),p_corrdiff(i)] = ttest2(hpc_pred{i}(~isnan(hpc_pred{i})), pfc_pred{i}(~isnan(pfc_pred{i})));
    [h_pred(i),p_pred(i)] = kstest2(hpc_pred{i}(~isnan(hpc_pred{i})), pfc_pred{i}(~isnan(pfc_pred{i})));
    xlim([0,0.75])
    
end
% 
%% 
%% Figure 5 - Showing the complexity not due to the other region being more complex
% 
% A: (difference in qOpt)

temp1 = zeros(3,nPartitions);
temp2 = zeros(3,nPartitions);
temp3 = zeros(3,nPartitions);

for j = 1:nPartitions
    for k = 1:nPatterns
        temp2(k,j) = Patterns(j,hpc,k).factorAnalysis.qOpt;
        temp3(k,j) = Patterns(j,pfc,k).factorAnalysis.qOpt;
        temp1(k,:) = temp2(k,:)./temp3(k,:);
    end
end
%
ratios = mean(temp1,2);
hpcQoptDims = mean(temp2,2);
pfcQoptDims = mean(temp3,2);


%% Figure 6
%% 
% *Pattern Specific* 

removePred
%  Remove tuples of (targetArea,pattern) from a given (targetArea, pattern)

figure;
g = gramm(...
    'x',rt.dimensionRemoved, ...
    'y', rt.performance, ...
    'color', categorical(rt.removePattern), ...
    'lightness', categorical(rt.sameDirectionLabel),...
    'linestyle', categorical(rt.sameDirectionLabel));
g = g.facet_grid(categorical(rt.basePatternLabel), categorical(rt.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar');
g = g.stat_summary('geom','area');
g.set_text_options("interpreter",'latex');
g = g.set_names('x','dims removed',...
    'y','peformance',...
    'column','Interaction',...
    'row','Pattern', ...
    'color','Pattern dims removed',...
    'lightness',"Remove Same/Different"+newline+"Pred. Target",...
    'linestyle',"Remove Same/Different"+ newline+"Pred. Target");
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
% Remove pattern from same target area

figure;
g = gramm(...
    'x',rt.dimensionRemoved, ...
    'y', rt.performance, ...
    'color', categorical(rt.removePattern), ...
    'subset',rt.sameDirection);
g = g.facet_grid(categorical(rt.basePatternLabel), categorical(rt.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar');
g = g.stat_summary('geom','area');
g.set_text_options("interpreter",'latex');
g=g.set_title("Removing dimensions" + newline + "(of similar target area only)");
g = g.set_names('x','dims removed',...
    'y','peformance',...
    'column','Interaction',...
    'row','Pattern', ...
    'color','Pattern dims removed');
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
%%
figure;
g = gramm(...
    'x',rt.dimensionRemoved, ...
    'y', rt.performance, ...
    'color', categorical(rt.removePattern), ...
    'subset',~rt.sameDirection);
g = g.facet_grid(categorical(rt.basePatternLabel), categorical(rt.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar');
g = g.stat_summary('geom','area');
g = g.set_color_options('lightness',200);
g.set_text_options("interpreter",'latex');
g = g.set_names('x','dims removed',...
    'y','peformance',...
    'column','Interaction',...
    'row','Pattern', ...
    'color','Pattern dims removed');
g=g.set_title('Removing different only')
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
% Overall target removal

figure;
g = gramm(...
    'x',rt.dimensionRemoved, ...
    'marker', rt.targetArea,...
    'y', rt.performance, ...
    'color', categorical(rt.sameDirectionLabel),'subset',~rt.sameDirection);
g = g.facet_grid([],categorical(rt.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar');
g = g.stat_summary('geom','area');
g = g.set_color_options('chroma',0,'lightness',30);
g.set_text_options("interpreter",'latex');
g = g.set_names('x','dims removed',...
    'y','peformance',...
    'column','Interaction',...
    'row','Pattern', ...
    'color','Brain area removed');
snapnow;
g=g.update('subset',rt.sameDirection);
%g = g.facet_grid([],categorical(rt.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar');
g = g.stat_summary('geom','area');
g = g.set_color_options(); % Restore default color
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
%% 
% Now let's see how similar those curves are in dimension reduced space

unstack(rt, '');

%% Figure 7
% 
% A: Dominant dimension predicitions

numUsedForPrediction = min(nTarget,nSource);
% make the averaged version
curr_cvLoss = cell(10,2,3);
curr_qOptDim = cell(10,2,3);
for p = 1:10
    for i = 1:nPatterns
        for j = 1:2
            if ~Patterns(p,j,i).singularWarning
                curr_cvLoss{p,j,i} = Patterns(p,j,i).factorAnalysis.cvLoss;
                curr_qOptDim{p,j,i} = Patterns(p,j,i).factorAnalysis.optDimFactorRegress;
            end
        end
    end
end
%%
figure(750)
clf
full_model_performance = [];

for i = 1:nPatterns
    for j = 1:2
        subplot(3,2,2*(i-1)+j)
        
        full_model = plots.plotPredictiveDimensions(numUsedForPrediction,curr_cvLoss(:,j,i),'optDim',curr_qOptDim, ...
            "mode", "fa");
        
        full_model_performance = [full_model_performance, full_model];
        xlim([0,12.5])
        hold on
        
        plot(1, full_model,'^');
        ylim([0, full_model+0.1])
        if j == 1
            ax1 = gca;
        else
            ax2 = gca;
        end
        title([Patterns(p,j,i).name Patterns(p,j,i).directionality])
    end
    linkaxes([ax1,ax2],'y')
    
end
