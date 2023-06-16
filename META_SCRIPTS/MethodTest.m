% Method test
addpath('/Volumes/MATLAB-Drive/');
addpath(genpath('/media/ryoung/Ark/SingleDayExpt'))
%                                                    
% Just a script for testing out a potential multilinear
% algebra approach to communication subpsaces in time
% (spaces that are allowed to change in amplitude over time
% and trials.)

%% Mode
waysOfPartitions = 2; %whether the script is bidirectional (split neurons in two ways)
run_selected_genH = false;
doFactorAnalysis  = false;
doTimeVarying = false;
doSingleNeuronPrediction = false;

%% Paths
if ~exist('pathSet','var') || pathSet == 0 % only set path if unset
    addpath(genpath(pwd)) % all folders in utils added, including semedo code
    
    % Whose computer are we running?
    if ispc
        paths = "C:\Users\BrainMaker\commsubspace\SingleDayExpt";
    elseif ismac || isunix
        paths = "/Volumes/sharespace-commsub/data";
    end
    % Set data paths for that computer
    arrayfun(@(path) addpath(genpath(path)), paths);
    % ---------- Add paths for specific users ------------------------------
    % addpath(genpath(datadef(user))); % this line throws an error for pc
    
    addpath(...
        fullfile(codedefine,"Shared"));
    addpath(genpath(fullfile(codedefine, "Shared", "utils")));
    pathSet = 1;
end

%% Script parameters
Default = struct();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Default.animal                        = "ZT2";
Default.generateH                     = "fromSpectra  "+" fromRipTimes";
Default.samplingRate                  = [] ;             % For spikes.getSpikeTrain, nan if not given
Default.winSize                       = {[-0.15, 0.15]}; % size of the window
Default.spikeBinSize                  = 0.030;           % 60 milliseconds, prviously 150 ms
Default.timesPerTrial                 = round(range(Default.winSize{1})/Default.spikeBinSize);
Default.sourceArea                    = "CA1";           % only when there are
Default.equalWindowsAcrossPatterns    = true;            % whether all three patterns have the same #windows
Default.singleControl                 = false;           % whether to use just one control column
Default.numPartition                  = 50;              % ways to split source and target
Default.binsToMatchFR                 = 20;
Default.quantileToMakeWindows         = 0.85;
Default.lowerControl                  = true;
Default.preProcess_FilterLowFR        = true;
Default.preProcess_matchingDiscreteFR = true;
Default.oldControlBehavior            = false;
Default.dimCompAnalysis               = 5;
Default.stablePerf                    = 0.9;
% Default.generateH = "fromWpli " + " fromRipTimes";
% Default.generateH = "fromCoherence "+" fromRipTimes";
Default.samplingRate  = [] ;              % For spikes.getSpikeTrain, nan if not given
Default.spikeBinSize  = 0.060;             % 60 milliseconds, prviously 150 ms
Default.timesPerTrial = 10;                % 10 times per trial
Default.winSize       = {[-0.15, 0.15]};     % size of the window
Default.sourceArea    = "CA1";             % only when there are
Default.equalWindowsAcrossPatterns = true;    % whether all three patterns have the same #windows
Default.singleControl = false;                 % whether to use just one control column
Default.numPartition = 5;                    % ways to split source and target
Default.binsToMatchFR = 20;
Default.quantileToMakeWindows = 0.7;
Default.lowerControl = true;
Default.preProcess_FilterLowFR = true;
Default.preProcess_matchingDiscreteFR = true;
Default.oldControlBehavior = false;

% Default.shortedRippleWindow = false; % using 100 milisecond window for ripples

% Write defaults to option struct if not provided
if ~exist('Option','var')
    Option = Default;
else
    for field = string(fieldnames(Default))'
        if ~isfield(Option, field)
            Option.(field) = Default.(field);
        end
    end
end

winSize = Option.winSize{1};
if contains(Option.generateH, "fromRipTimes")
    Option.generateFromRipTimes = true;
else
    Option.generateFromRipTimes = false;
end

animal = Option.animal;
disp("Running with Option struct => ")
disp(Option);
%% Shortcut/alias variables to improve readability
THETA  = 1;
DELTA  = 2;
RIPPLE = 3;
if Option.sourceArea == "CA1"
    HPC = 1;
    PFC = 2;
else
    PFC = 1;
    HPC = 2;
end
patternNames = ["theta","delta","ripple"];

%% Mung/Clean the data
frequenciesPerPattern = [6 10; 0.5 4; 150 200];
[nPatterns,~] = size(frequenciesPerPattern);

if Option.singleControl == true
    nPatternAndControl = nPatterns+1;
else
    nPatternAndControl = nPatterns*2;
end


%%%%%%%%%%%%%%%% READ IN RAW PATTERN DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find network pattern events
if contains(Option.generateH, "fromSpectra")
    load(Option.animal + "spectralBehavior.mat");
    if Option.sourceArea == "CA1"
        spectrogram   = efizz.S1;
    else
        spectrogram   = efizz.S2;
    end
    
    frequencyAxis = efizz.f;
    Htimes = efizz.t;
    [H, Hvals, Hnanlocs, Htimes] = eventMatrix.generateFromSpectra(Htimes, spectrogram, frequencyAxis,...
        frequenciesPerPattern);
elseif contains(Option.generateH, "fromFilteredEEG")
    load(Option.animal + "avgeeg.mat");
    [runningSessions, sleepSessions] = getRunningSessions(Option.animal);
    [H, Hvals, Hnanlocs, Htimes] = eventMatrix.generateFromFilteredEEG(avgeeg, ...
        Option.sourceArea, "patterns",patternNames(1:3),"downsample",10, "sleepSessions", sleepSessions);
elseif contains(Option.generateH, "fromCoherence")
    load(Option.animal + "spectralBehavior.mat");
    spectrogram = efizz.C;
    frequencyAxis = efizz.f;
    Htimes = efizz.t;
    [H, Hvals, Hnanlocs, Htimes] = eventMatrix.generateFromSpectra(Htimes, spectrogram, frequencyAxis,...
        frequenciesPerPattern);
elseif contains(Option.generateH, "fromWpli")
    load(Option.animal + "spectralBehavior.mat");
    spectrogram = efizz.wpli;
    frequencyAxis = efizz.f;
    Htimes = efizz.t;
    [H, Hvals, Hnanlocs, Htimes] = eventMatrix.generateFromSpectra(Htimes, spectrogram, frequencyAxis,...
        frequenciesPerPattern);
    
else
    error("Core method for deriving the event matrix is not recognized");
end

%% Modify ripple pattern? Lower the threshold as with window sizes
if contains(Option.generateH,"fromRipTimes")
    load(Option.animal + "globalripple01.mat");
    
    if any(contains(Option.generateH, ["fromWpli", "fromCoherence"]))
        [~, H(:,RIPPLE), Hnanlocs(:,RIPPLE), Hvals(:,RIPPLE), minRippleThreshold, original] = ...
            eventMatrix.generateFromRipples(globalripple, ...
            'amplitude_at_riptime', true,...
            'rippleBand', Hvals(:,RIPPLE),... RY: Hvals, not H here for obvious reasons: you  want the original ripple band activity
            'rippleBandTime', Htimes,...
            'globalrippleWindowUnits', 'amp');
    else
        [~, H(:,RIPPLE), Hnanlocs(:,RIPPLE), Hvals(:,RIPPLE), minRippleThreshold, original] = ...
            eventMatrix.generateFromRipples(globalripple, ...
            'amplitude_at_riptime', true,...
            'rippleBand', Hvals(:,RIPPLE),... RY: Hvals, not H here for obvious reasons: you  want the original ripple band activity
            'rippleBandTime', Htimes,...
            'globalrippleWindowUnits', 'std');
    end
end

%%%%%%%%%%%%%%%% WINDOW SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Making Windows
% windows of network patterns
[cellOfWindows, cutoffs] = windows.make(   Htimes, Option.quantileToMakeWindows, H(:,THETA:DELTA), winSize);
if any(contains(Option.generateH, ["fromCoherence","fromWpli"]))
    [cellOfWindows(RIPPLE), cutoffs(RIPPLE)] = windows.make(Htimes,...
        Option.quantileToMakeWindows,  H(:,RIPPLE), winSize,...
        'quantile', Hvals(:,RIPPLE),'higherThanQuantile', true); % % RY: quantile needs to be hvals for ripple coherence/wpli threshold to be correct, but timesd computed from H such that non-ripple times thrown out
else
    [cellOfWindows(RIPPLE), cutoffs(RIPPLE)] = windows.make(Htimes, ...
        1,   H(:,RIPPLE), winSize, 'threshold', 'raw','higherThanQuantile', true);
end
    
% end
disp(newline);
disp('--------------------------')
disp('Initiale window  creation:')
disp('--------------------------')
cellfun(@(x,y)  fprintf("%d windows for %s\n", size(x,1),  y), ...
    cellOfWindows, cellstr(patternNames(1:3)));
cellfun(@(x,y)  fprintf("%d timerange for %s\n", range(x(:)),  y), ...
    cellOfWindows, cellstr(patternNames(1:3)));

% equalize number of windows across patterns based on input argument
if Option.equalWindowsAcrossPatterns == true
    cellOfWindows = windows.equalizeWindowsAcrossPatterns(cellOfWindows, nPatterns);
end

disp(newline);
disp('--------------------------')
disp('Equalized windows        :')
disp('--------------------------')
cellfun(@(x,y)  fprintf("%d windows for %s\n", size(x,1),  y), ...
    cellOfWindows, cellstr(patternNames(1:3)));

numWindowsCut = size(cellOfWindows{1},1);
%windows.printWindowOverlap(cellOfWindows, patternNames);

%%%%%%%%%%%%%%%% CONTROL SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All non-ripple patterns
% -----------------------
if Option.lowerControl
    quantileControl = 1 - Option.quantileToMakeWindows; %TODO at some point we may want an Option.(something) to toggle this
else
    quantileControl = Option.quantileToMakeWindows;
end

if Option.oldControlBehavior
    Hc =  control.generatePatternShuffle(H(:,1:3), Htimes, cellOfWindows); % add control patterns;
else
    Hc = H;
end
[Hc_cellOfWindows, Hc_cutoffs] = windows.make(Htimes,  quantileControl,...  % add windows of control patterns
    Hc(:,THETA:DELTA), winSize,... % Selects less than quantile
    'higherThanQuantile', Option.oldControlBehavior);

% -------
% Ripples
% -------
if any(contains(Option.generateH, ["fromCoherence","fromWpli"]))
    [Hc_cellOfWindows(RIPPLE), Hc_cutoffs(RIPPLE)] = windows.make(Htimes,...
        quantileControl,  Hc(:,RIPPLE), winSize,...
        'quantile', Hvals(:,RIPPLE),'higherThanQuantile', Option.oldControlBehavior);
else
    if Option.oldControlBehavior
        [Hc_cellOfWindows(RIPPLE), Hc_cutoffs(RIPPLE)] = windows.make(Htimes, ...
            1,   Hc(:,RIPPLE), winSize, 'threshold', 'raw','higherThanQuantile', true); %RY quantile won't work because these are raw
    else
        [Hc_cellOfWindows(RIPPLE), Hc_cutoffs(RIPPLE)] = windows.make(Htimes, ...
            quantileControl,   Hvals(:,RIPPLE), winSize, 'threshold', 'quantile','higherThanQuantile', Option.oldControlBehavior); %RY quantile won't work because these are raw
    end
end

% -----
% Clean
% -----
% clean up control windows: remove each control pattern's window's overlap
for pattern = 1:nPatterns
    curr = windows.removeOverlapsBetweenPattern(...
        cell2mat(cellOfWindows(:,pattern)), cell2mat(Hc_cellOfWindows(:,pattern)));
    Hc_cellOfWindows{pattern} = curr;
end

% -----------------------------------
% Merge many controls into 1 control?
% -----------------------------------
% % Merge into one
cellOfWindows(nPatterns+1:nPatterns*2) = Hc_cellOfWindows;
cutoffs = [cutoffs,Hc_cutoffs];

% -----------------------------------------
% Ensure each pattern has equal # of window
% -----------------------------------------
% Equalize trials/windows for each pair of patttern-controlPattern
[cellOfWindows, warnedEmptyControls] =...
    control.equalizePatternControl(nPatterns, cellOfWindows);

% pick control pattern that actually contains controls, would break if all
% three are empty...
if Option.singleControl
    if warnedEmptyControls
        for iPossibleControl = nPatterns+1:nPatterns*2
            if ~isempty(cellOfWindows{iPossibleControl})
                cellOfWindows{nPatterns+1} = cellOfWindows{iPossibleControl};
                disp("here")
            end
        end
    end
end

if ~any(contains(patternNames,"control"))
    patternNames = [patternNames; patternNames+"-control"]';
    patternNames = patternNames(:)';
end

%%%%%%%%%%%%%%%% SPIKE SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Getting spikes
[timeBinStartEnd, timeBinMidPoints, ~, spikeCountMatrix, spikeRateMatrix, ...
    areaPerNeuron, cell_index, sessionTypePerBin] = spikes.getSpikeTrain(Option.animal, ...
    Option.spikeBinSize,  Option.samplingRate);

if Option.preProcess_FilterLowFR % filter the neurons whose firing rate is lower than specified threshold
    [spikeCountMatrix, spikeRateMatrix, avgFR, areaPerNeuron, cell_index]...
        = trialSpikes.filterFR(spikeCountMatrix, spikeRateMatrix, 0.1, ...
                               timeBinStartEnd, areaPerNeuron, cell_index);
end

[spikeSampleMatrix, spikeSampleTensor, trialTimes] = trialSpikes.generate(...
    spikeCountMatrix,...
    timeBinStartEnd, cellOfWindows, ...
    Option.timesPerTrial, nPatternAndControl);
%% Separate spikesSampleMatrix/Tensor by area that neurons are in PFC and neurons that in HPC
[pfcFR, hpcFR] = trialSpikes.separateFiringRate(avgFR, areaPerNeuron);
T_pfc = trialSpikes.separateSpikes(spikeSampleTensor, areaPerNeuron, "PFC");
T_hpc = trialSpikes.separateSpikes(spikeSampleTensor, areaPerNeuron, "CA1");
%% Separate firing pattern into source and target
[nPFCneurons,~,~] = size(T_pfc{1});
[nHPCneurons,~,~] = size(T_hpc{1});

%% Create the tensor subspace of firing outer products
celllookup = cellInfo.getCellIdentities(animal, cell_index, areaPerNeuron);


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
    [beh,badtimes] =  table.behavior.lookup(animal, timeBinMidPoints, []);
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
