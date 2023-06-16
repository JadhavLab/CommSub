clear all

% ---- PATH -----
% Matlab uses startup.m to run startup code...
% Put this in your startup.m so that the code for this is in path:
%
% addpath(genpath('/Volumes/MATLAB-Drive/')) % or wherever your CODE files are located
% addpath(genpath('~/Data/Raw/')) % or wherever your DATA files are located

% ===================================================
% OPTION STRUCT encoding properties of the script run
% ===================================================
% see +option.default() to set default options
if ~exist('Option','var')
    Option = option.defaults(); 
else
    Option = option.setdefaults(Option);
end

%%%%%%%%%%%%%%%% DISPLAY OUR OPTIONS TO USER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Running with Option struct => ")
disp(Option);

%%%%%%%%%%%%%%%% OBTAIN EVENT MATRICES    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Events] = events.ThetaDeltaRipple(Option);

%%%%%%%%%%%%%%%% CUT WINDOWS WITH EVENT MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[cellOfWindows, cutoffs] = windows.ThetaDeltaRipple(Events, Option);

%%%%%%%%%%%%%%%% ACQUIRE SPIKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Getting spikes
[timeBinStartEnd, timeBinMidPoints, ~, spikeCountMatrix, spikeRateMatrix, ...
    areaPerNeuron, cell_index, sessionTypePerBin] = spikes.getSpikeTrain(Option.animal, ...
    Option.spikeBinSize,  Option.samplingRate);

% filter the neurons whose firing rate is lower than specified threshold
if Option.preProcess_FilterLowFR 
    [spikeCountMatrix, spikeRateMatrix, avgFR, areaPerNeuron, cell_index]...
        = trialSpikes.filterFR(spikeCountMatrix, spikeRateMatrix, 0.1, ...
        timeBinStartEnd, areaPerNeuron, cell_index);
end

%%%%%%%%%%%%%%%% ACQUIRE TRIALS FROM WINDOWS + SPIKES %%%%%%%%%%%%%%%%%%%
[spikeSampleMatrix, spikeSampleTensor, trialTimes] = trialSpikes.generate(...
    spikeCountMatrix, timeBinMidPoints, cellOfWindows, ... % RYAN bug here .. timeBinStartEnd instead of timeBinMidPoints
    Option.timesPerTrial, Option.nPatternAndControl);


% %%%%%%%%%%%%%%% SETUP RAW DATA STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure for separated data
r.hpc = struct;
r.pfc = struct;
r.avgFR = avgFR;
r.areaPerNeuron = areaPerNeuron;
r.windowInfo.cellOfWindows = cellOfWindows;
r.windowInfo.nWindows      = cellfun(@(x) size(x, 1), cellOfWindows);
r.nPattern = Option.nPattern;
r.nControl = nPatternAndControl - r.nPattern;

%%%%%%%%%%%%%%%% SEPRATE BRAIN AREA DATASULT STRUCTURES %%%%%%%%%%%%%%%%%%
%% Separate spikesSampleMatrix/Tensor by area that neurons are in PFC and neurons that in HPC
[r.pfc.FR, r.hpc.FR] = trialSpikes.separateFiringRate(avgFR, areaPerNeuron);
r.pfc.T = trialSpikes.separateSpikes(spikeSampleTensor, areaPerNeuron, "PFC");
r.hpc.T = trialSpikes.separateSpikes(spikeSampleTensor, areaPerNeuron, "CA1");
r.pfc.X = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "PFC");
r.hpc.X = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "CA1");

%% Separate firing pattern into source and target
[nPFCneurons,~,~] = size(r.pfc.X{1});
[nHPCneurons,~,~] = size(r.hpc.X{1});
r.pfc.nNeurons = nPFCneurons;
r.hpc.nNeurons = nHPCneurons;

r.celllookup = cellInfo.getCellIdentities(Option.animal, cell_index,...
                                          areaPerNeuron);

%%%%%%%%%%%%%%%% SETUP PARTITIONS AND RESULT STRUCTURES %%%%%%%%%%%%%%%%%%
Patterns = trialSpikes.partitionAndInitialize(r, Option);


% Variables to produce
% spikeRateMatrix, timeBinMidPoints, 
% Hvals, Htimes, efizz, avgeeg

%% -----------------
%% TENSOR PREPROCESS
%% -----------------
%method = ["spikerate-areawiseOuterProduct-flatten2-compress", "", ""]; % Look for patterns in triple cumulant of spikes_area1, spikes_area2, lfp
method = ["spikerate", "", ""]; % Look for patterns in triple cumulant of spikes_area1, spikes_area2, lfp
decomposition_method = "preferential-subspace-identification";
hankelize_time = 1; % 1 = no hankelization
rank_setting = 0;

%method = ["spectral", "spikerate", "spikes*lfp"]; % Look for patterns in triple cumulant of spikes_area1, spikes_area2, lfp

% Tensor commsubspace can get so large that it has to be majorly down-sampled even fit.

% Tensor of spikes?
% -----------------
spike_method = method(1);
if contains(spike_method, "spikerate")
    % Compute tensor for all data, annotated with lfp
    disp("Spikerate")
    [T, T_times] = tensor.spikes(spikeRateMatrix, timeBinMidPoints, hankelize_time,  'Dim', 2);
    if contains(spike_method, "areawiseOuterProduct")
        disp("->outerProduct")
        T = tensor.areawiseOuterProduct(T, areaPerNeuron, 'big', true);
    end
    if contains(spike_method, "diagSqr")
        disp("->diagSqr")
        T = permute(T,[3,1,2]);
        c = size(T, 1);
        idx = 1:c+1:numel(T(1,:,:));
        for t = 1:size(T,3)
            T(t,idx) = sqrt(T(t, idx));
        end
        T = ipermute(T,[3,1,2]);
    end
    keyboard
    if contains(spike_method, "flatten2")
        disp("->flatten")
        T = permute(T, [3, 1 2]);
        T = reshape(T, size(T,1), []);
        T = T';
    end
    if contains(spike_method, "compress")
        n = 300;
        T = T';
        mu = mean(T,1);
        [score, coeff, latent] = pca(T-mu);
        T = coeff(:, 1:n);
    end
else
    T = []
end

% Tensor of LFP?
% --------------
lfp_method   = method(2);
if lfp_method == "H"
    [Th, Th_times] = tensor.H(Hvals, Htimes, hankelize_time);
elseif lfp_method == "spectral"
    [Th, Th_times] = tensor.spectral(efizz, hankelize_time, 'fields', ["S1","S2","wpli"]);
elseif lfp_method == "hilbert"
    [Th, Th_times] = tensor.hilbert(avgeeg, hankelize_time);
else
    Th = [];
end

% Final tensor form
% -----------------
final_method = method(3);

% ... same times?
if all(contains(final_method, ["lfp","spikerate"]))
    [T, Th, times] = tensor.unifyTimes(T, T_times, Th, Th_times);
end

% ... apply combine method
if final_method == "lfp+spikerate" % combine spike and lfp in a single tensor mode
    % LFP AND SPIKES SEPARATE OUTER PHENOMENA
    T = cat(1, T, Th);
elseif final_method == "lfp*spikerate"% product, place spike and lfp in different tensor modes
    % LFP AND SPIKES COMBINED INNER/CUMULANT PHENOMENA
    [T, times] = tensor.outerProduct(T, Th, times, 'downsample', -1, 'center', false);
elseif final_method == "lfp"
    % JUST LFP
    T = Th;
    times = Th_times;
elseif final_method == "spikes"
    % JUST SPIKES
    % T = T;
    times = T_times;
elseif final_method == "use_existing_windows" % use subspace windows of activity for final tensor
    % SEPARATE SETS OF SPIKES PER EVENT
    [T] = tensor.areawiseOuterProduct(T_hpc, T_pfc);
end

% --------------------------------------------------------------
% Now that we have our tensors, we
% 1. find the rank
% 2. apply a type of tensor decomposition to learn the factors
% --------------------------------------------------------------

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%% CANNONICAL TENSOR DECOMPOSITION
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
if decomposition_method == "cpd"
    if rank_setting == 0
        % Compute rank
        RankEst = repmat(struct(),6,1);
        parfor i = 1:numel(T)
            [RankEst(i).R, RankEst(i).L_lb, RankEst(i).L_cpd] = rankest(T{i});
        end
    end

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%% MULTILINEAR TENSOR DECOMPOSITION
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
% Obtain the multilinear tensor rank utilizing tensorlab
% ------------------------------------------------------
elseif decomposition_method == "mlsvd"

    if rank_setting == 0
        % (Takes a looooong time)
        % (https://www.tensorlab.net/doc/)
        MlRankEst = repmat(struct(),6,1);
        parfor i = 1:numel(T)
            [MlRankEst(i).R, MlRankEst(i).L_lb, MlRankEst(i).L_cpd] = mlrankest(T{i});
        end
    end

elseif decomposition_method == "data_fusion_method"

    if combine_method ~= "use_existing_windows"
        error("Must use use_existing_windows for datafusion");
    end

    % Get rank?
    if rank_setting == 0
        error("Not implemented. Please provide rank.");
    end
    

    % Compute UV matrix factors
    allPatternsSubspaceSimilarity; 
    dim = min(size(B_{1})); % chose the top five dimensions
    celllookup = cellInfo.getCellIdentities(animal, cell_index, areaPerNeuron);
    [subspace_belonging, cell_subspaces] = components.subspaceSimilarity(dim, B_, spikeRateMatrix, celllookup,"CA1");

    % Compute tensor (UV) factors
    for i = 1:numel(T)
        [decomp{i}, model{i}] = components.UVtensorDecomp(T{i}, B_{i}, optDim{i}, rank_setting)
    end

    % Visualize factors found
    fig('Factor visualization');
    factors = string(fieldnames(sol.factors));
    tiledlayout(10, numel(factors), 'TileSpacing', 'compact', 'Padding', 'compact'); 
    colormaps = {'buda', 'acton', 'batlow', 'nuuk','tokyo'};
    for i = 1:min(size(sol.factors.(factors(f)),2),10)
        for f = 1:numel(factors)
            colors = crameri(colormaps{f}, min(size(sol.factors.T,2),10)); 
            nexttile;
            if factors(f).lower() == "t"
                X = smoothdata(sol.factors.(factors(f))(:,i), 'SmoothingFactor',0.80);
                kws = {'linestyle','-'};
            else
                X = sol.factors.(factors(f))(:,i);
                kws = {'marker','o', 'linestyle', 'none','markerfacecolor','auto', 'markersize',4};
            end
                plot(X,'Color',colors(i,:), kws{:});
            hold on;
            set(gca,'xtick',[]);
        end
    end

     tensor_behavior = arrayfun(@(x) table.behavior.lookup(animal, mean(cellOfWindows{x},2), [], 'valueOnly', true, 'scaleVars', true), 1:numel(patternNames), 'UniformOutput', false);
     [matrix_behavior, badtimes] = table.behavior.lookup(animal, timeBinMidPoints, [], 'valueOnly', true, 'scaleVars', true);
     uvTimes = timeBinMidPoints(~badtimes);
     cell_subspaces = cellfun(@(x) x(:,~logical(badtimes))', cell_subspaces, 'UniformOutput', false);
     
     f=fig('Tensor Behavior Relevance'); clf;
     tiledlayout(f,2,3);
     arrayfun(@(i)tensor.findTimeFactorBehaviorCorr(decomp{i}, 'T', tensor_behavior{i},'colormap',colormaps{i},'constrainSig',0.01), 1:3, 'UniformOutput', false);
     %f=gcf;
     %for i = 1:3
     %    title(f.Children.Children(3-i + 1), patternNames(i));
     %end
     sgtitle("Network-Pattern-Specific" + newline + "Communication Subspace Modes" + newline + "Correlations WITH BEH" + newline + "(Only p < 0.01)");
     %fig("Matrix Behavior Relevance"); 
     arrayfun(@(i) tensor.findTimeFactorBehaviorCorr(cell_subspaces{i}, [], matrix_behavior,'colormap',colormaps{i},'constrainSig',0.01), 1:3, 'UniformOutput', false);
     %f=gcf;
     %for i = 1:3
     %    title(f.Children.Children(3-i + 1), patternNames(i));
     %end
     %sgtitle("Network-Pattern-Specific" + newline + "Communication Subspace Modes" + newline + "Correlations WITH BEH" + newline + "(Only p < 0.01)");

elseif contains(decomposition_method, "preferential-subspace-identification")

    % Data
    T = T';
    y = T;
    [beh,badtimes] =  table.behavior.lookup(Option.animal, timeBinMidPoints, []);
    y(badtimes,:) = [];

    % Settings
    beh_fields = ["rewarded","trajbound","directional_lindist","leftright", "X","Y"];
    nx = 6;
    ni = min(3, numel(beh_fields));
    if contains(spike_method,'outer')
        hankel_embedding = 8;
    else
        hankel_embedding = 20; % the hankelization horizon (5 time points delay embedding), for creating block hankel matrix
    end

    % Calc psid on different sets
    psid_analysis = struct();
    params = {beh, nx, ni, hankel_embedding, 'beh_fields', beh_fields, 'train', "neighbor"};
    disp("Calculating total")
    [psid_analysis.total, z, B] = psid.calculate(y, params{:}, 'set_name', 'total', 'bigdata', true);

    params = {beh, 4, ni, hankel_embedding, 'beh_fields', beh_fields, 'train', "neighbor"};
    disp("Calculating CA1")
    psid_analysis.ca1   = psid.calculate(y(:,areaPerNeuron == "CA1"), params{:}, 'set_name', 'ca1', 'bigdata', true);
    disp("Calculating PFC")
    psid_analysis.pfc   = psid.calculate(y(:,areaPerNeuron == "PFC"), params{:}, 'set_name', 'pfc', 'bigdata', true);

    %% Plot results : Corr
    %sets = ["ca1","pfc"];
    sets = ["total"];
    for psid_set = sets
        psid.plots(psid_analysis.(psid_set), beh, 'beh_fields', beh_fields, 'overall_projLowerD', false);
    end

    %% Annotate behavior with space data
    B=[];
    for field = sets
        B.(field) = psid.annotate(psid_analysis.(field), beh, 'tag',field,'as','rows');
    end

    %% Avg traj
    traj = [];
    for field = sets
        fig(field + " traj average");
        traj.(field) = psid.trajavg(B.(field), 80, 'ploton', true);
        sgtitle("PSID " + field);
    end

   %%  
   mkdir('~/Data/psid')
   hashkey = DataHash([animal, method, decomposition_method]);
   hashkey = hashkey(1:8); 
   mfile = "~/Data/psid/psid_analysis_" + hashkey;
   dfile = "~/Data/psid/psid_data_" + hashkey;
   save(mfile, '-struct', 'psid_analysis','-v7.3');
   %save(dfile, 'y', 'beh', '-v7.3');
   writetable(beh, '~/Data/psid/psid_behavior.csv')
   vars = num2cell([animal, method, decomposition_method, hashkey]); 
   varnames = ["animal", "spike", "lfp", "spikelfpcomb", "decomposition", "hash"];
   tab = table(vars{:}, 'VariableNames', varnames);
   writetable(tab, '~/Data/psid/psid_analysis_table.csv', 'WriteMode', 'append');

end
