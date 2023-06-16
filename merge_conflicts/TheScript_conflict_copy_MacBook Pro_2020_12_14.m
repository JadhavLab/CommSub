
%% Mode
waysOfPartitions = 2; %whether the script is bidirectional (split neurons in two ways)

run_selected_genH = false;
%% Paths
if ~exist('pathSet','var') || pathSet == 0 % only set path if unset
    addpath(genpath(pwd)) % all folders in utils added, including semedo code
    
    % Whose computer are we running?
    if ispc
        paths = " C:\Users\BrainMaker\commsubspace\SingleDayExpt";
    elseif ismac
        paths(1) = "/Volumes/sharespace-commsub/data";
        paths(2) = "~/Data/commsubspace";
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
% determined by spectral power with pattern determined by global ripple
Default.animal    = "ZT2";
% Default.generateH = "fromSpectra "+" fromRipTimes";
Default.generateH = "fromFilteredEEG "+" fromRipTimes";
% Default.generateH = "fromWpli " + " fromRipTimes";
% Default.generateH = "fromCoherence "+" fromRipTimes";
Default.samplingRate  = [] ;              % For spikes.getSpikeTrain, nan if not given
Default.spikeBinSize  = 0.15;               % 100 milliseconds
Default.timesPerTrial = 10;                % 10 times per trial
Default.winSize       = {[-0.15, 0.15]};     % size of the window
Default.sourceArea    = "CA1";             % only when there are
Default.equalWindowsAcrossPatterns = true;    % whether all three patterns have the same #windows
Default.singleControl = false;                 % whether to use just one control column
Default.numPartition = 50;                    % ways to split source and target
Default.binsToMatchFR = 20;
Default.quantileToMakeWindows = 0.7;
Default.lowerControl = true;
Default.preProcess_FilterLowFR = true;
Default.preProcess_matchingDiscreteFR = true;
Default.shortedRippleWindow = false; % using 100 milisecond window for ripples

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

%%%%%%%%%%%%%%%% READ IN RAW PATTERN DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
animal = Option.animal;
disp("Running with Option struct => ")
disp(Option);
%% Shortcut/alias variables to improve readability
THETA = 1;
DELTA = 2;
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
    [~, sleepSessions] = getRunningSessions(Option.animal);
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
if Option.shortedRippleWindow
    rippleWinSize = [-0.05, 0.05];
    if any(contains(Option.generateH, ["fromCoherence","fromWpli"]))
        [cellOfWindows(RIPPLE), cutoffs(RIPPLE)] = windows.make(Htimes,...
            Option.quantileToMakeWindows,  H(:,RIPPLE),  rippleWinSize,...
            'quantile', Hvals(:,RIPPLE), 'higherThanQuantile', true);
    else
        [cellOfWindows(RIPPLE), cutoffs(RIPPLE)] = windows.make(Htimes, ...
            1,   H(:,RIPPLE), rippleWinSize, 'threshold', 'raw','higherThanQuantile', true);
    end
else
    if any(contains(Option.generateH, ["fromCoherence","fromWpli"]))
        [cellOfWindows(RIPPLE), cutoffs(RIPPLE)] = windows.make(Htimes,...
            Option.quantileToMakeWindows,  H(:,RIPPLE), winSize,...
            'quantile', Hvals(:,RIPPLE),'higherThanQuantile', true); % % RY: quantile needs to be hvals for ripple coherence/wpli threshold to be correct, but timesd computed from H such that non-ripple times thrown out
    else
        [cellOfWindows(RIPPLE), cutoffs(RIPPLE)] = windows.make(Htimes, ...
            1,   H(:,RIPPLE), winSize, 'threshold', 'raw','higherThanQuantile', true);
    end
    
end
disp(newline);
disp('Initiale window  creation:')
disp('--------------------------')
cellfun(@(x,y)  fprintf("%d windows for %s\n", size(x,1),  y), ...
    cellOfWindows, cellstr(patternNames(1:3)));

% equalize number of windows across patterns based on input argument
if Option.equalWindowsAcrossPatterns == true
    cellOfWindows = windows.equalizeWindowsAcrossPatterns(cellOfWindows, nPatterns);
end

disp(newline);
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

% Clean
% -----
% clean up control windows: remove each control pattern's window's overlap
for pattern = 1:nPatterns
    curr = windows.removeOverlapsBetweenPattern(...
        cell2mat(cellOfWindows(:,pattern)), cell2mat(Hc_cellOfWindows(:,pattern)));
    Hc_cellOfWindows{pattern} = curr;
end

% Merge many controls into 1 control?
% -----------------------------------
% % Merge into one
cellOfWindows(nPatterns+1:nPatterns*2) = Hc_cellOfWindows;
cutoffs = [cutoffs,Hc_cutoffs];

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
[timeAxis, times_spiking, spikeCountMatrix, spikeRateMatrix, ...
    areaPerNeuron, cell_index] = spikes.getSpikeTrain(Option.animal, ...
    Option.spikeBinSize,  Option.samplingRate);

if Option.preProcess_FilterLowFR % filter the neurons whose firing rate is lower than specified threshold
    [spikeCountMatrix, spikeRateMatrix, avgFR, areaPerNeuron, cell_index]...
        = trialSpikes.filterFR(spikeCountMatrix, spikeRateMatrix, 0.1, ...
        timeAxis, areaPerNeuron, cell_index);
end

[spikeSampleMatrix, spikeSampleTensor] = trialSpikes.generate(spikeCountMatrix,...
    timeAxis, cellOfWindows, ...
    Option.timesPerTrial, nPatternAndControl);
%% Separate spikesSampleMatrix/Tensor by area that neurons are in PFC and neurons that in HPC
[pfcFR, hpcFR] = trialSpikes.separateFiringRate(avgFR, areaPerNeuron);
X_pfc = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "PFC");
X_hpc = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "CA1");

%% Separate firing pattern into source and target
[nPFCneurons,~] = size(X_pfc{1});
[nHPCneurons,~] = size(X_hpc{1});

%%%%%%%%%%%%%%%% PRE-RESULT SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results place to store outputs
clear Patterns
Patterns = struct("X_source",[], "X_target",[]);
Patterns.rankRegress = struct(...
    "B", [], ...
    "B_", [], ...
    "optDimReducedRankRegress", 0, ...
    "singlesource_B", [], ...
    "singlesource_optDim",[]);
Patterns.factorAnalysis = struct("qOpt", []);

%---- experimental code starts here
patternNames = ["theta","delta","ripple",...
    "theta-control","delta-control","ripple-control"];
if nPatternAndControl == nPatterns+1
    patternNames = ["theta","delta","ripple","control"];
end
Patterns = repmat(Patterns, ...
    [Option.numPartition, waysOfPartitions, nPatternAndControl]);


for iPartition = 1:Option.numPartition
    
    % --------------------------------------------------------------------------
    % Execute partitioning or firing rate matching
    % --------------------------------------------------------------------------
    % TODO COMPRESS THESE INTO SINGLE METHOD?
    switch waysOfPartitions
        case 4 % 4 WAY SPLIT
            [s_hpc, s_pfc, t_hpc, t_pfc, s_hpc_index, s_pfc_index,...
                t_hpc_index,t_pfc_index ] = trialSpikes.twoWaySplit...
                (pfcFR, hpcFR, X_hpc, X_pfc, Option.binsToMatchFR);
            directionality = ["hpc-hpc","hpc-pfc", "hpc-hpc", "pfc-hpc"];
        case 3 % 3 WAY SPLIT (three-way, where source split in two)
            if Option.targetArea == "CA1"
                [s_hpc, s_pfc, t_ca1, s_hpc_index, s_pfc_index, t_ca1_index] ...
                    = trialSpikes.threeGroupsSplit(hpcFR, pfcFR, X_hpc, X_pfc, Option.binsToMatchFR);
            else
                [s_pfc, s_hpc, t_pfc, s_pfc_index, s_hpc_index, t_pfc_index] ...
                    = trialSpikes.threeGroupsSplit(pfcFR, hpcFR, X_pfc, X_hpc, Option.binsToMatchFR);
            end
            if targetArea == "CA1"
                directionality = ["hpc-hpc","pfc-hpc"];
            else
                directionality = ["pfc-hpc", "pfc-pfc"];
            end
        case 2 % TWO WAY SPLIT (three-way, where target split in two)
            if Option.sourceArea == "CA1"
                directionality = ["hpc-hpc","hpc-pfc"];
            else
                directionality = ["pfc-hpc","pfc-pfc"];
            end
            if ~Option.preProcess_matchingDiscreteFR
                % (DEPRICATED) : not really in use as much
                % randomly select neurons based on number of neurons in the
                % respective regions and split source and target
                [X_source,X_target, nSource, nTarget,...
                    index_source, index_target] = ...
                    trialSpikes.splitSourceTarget...
                    (2, nPatternAndControl, X_hpc, X_pfc, 'specifiedSource', Option.sourceArea);
            else
                % (PRIMARY METHOD)
                % match the firing rate distribution
                if lower(Option.sourceArea) == "ca1"
                    [s_hpc, target, nSource, nTarget,...
                        s_hpc_index, index_target] = ...
                        trialSpikes.matchFRinDiscreteRanges...
                        (X_hpc, X_pfc, hpcFR, pfcFR, Option.binsToMatchFR);
                    t_hpc = target(1,:);
                    t_pfc = target(2,:);
                    t_hpc_index = index_target(1,:);
                    t_pfc_index = index_target(2,:);
                elseif lower(Option.sourceArea) == "pfc"
                    [s_pfc,target, nSource, nTarget,...
                        s_pfc, index_target] = ...
                        trialSpikes.matchFRinDiscreteRanges...
                        (X_hpc, X_pfc, hpcFR, pfcFR, Option.binsToMatchFR);
                    t_pfc = target(1,:);
                    t_hpc = target(2,:);
                    t_pfc = index_target(1,:);
                    t_hpc = index_target(2,:);
                end
            end
        otherwise
            error('Option not implemented')
    end
    % --------------------------------------------------------------------------
    
    % ------------------------------
    % Place paritioned data properly
    % ------------------------------
    for i = 1:numel(patternNames)
        for j = 1:numel(directionality)
            
            % Parse directionality
            sourcetarg = directionality.split('-');
            source = sourcetarg(1);
            target = sourcetarg(1);
            if source == "hpc"
                s_dat = s_hpc;
                s_ind = s_hpc_index;
            elseif source == "pfc"
                s_dat = s_pfc;
                s_ind = s_pfc_index;
            end
            if target == "hpc"
                t_dat = t_hpc;
                t_ind = t_hpc_index;
            elseif target == "pfc"
                t_dat = t_pfc;
                t_ind = t_pfc_index;
            end
            
            % Assign x_source and x_target
            Patterns(iPartition,j,i).X_source = s_dat{i};
            Patterns(iPartition,j,i).X_target = t_dat{i};
            
            % Assign index_source and index_target
            Patterns(iPartition,j,i).index_source = s_ind;
            Patterns(iPartition,j,i).index_target = t_ind;
            
            % Assign directionality
            Patterns(iPartition,j,i).directionality = directionality(j);
            
            % Assign pattern name
            Patterns(iPartition,j,i).name = patternNames(i);
        end
    end
    
end
% experimeental code stops here
%%
%%%%%%%%%%%%%%%% RANK-REGRESS SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if waysOfPartitions ~= 2
    nTarget = size(Patterns(1,1,1).X_target,1);
    nSource = min(size(Patterns(1,1,1).X_source,1),...
        size(Patterns(1,2,1).X_source,1));
end
numDimsUsedForPrediction = 1:min(nTarget,nSource);

for p = 1:Option.numPartition
    disp("processing rrr for "+p+" partition")
    for i = 1:nPatternAndControl
        
        B_singleprediction = cell(1,nSource);
        dim_singleprediction = cell(1,nSource);
        
        % Number of cross validation folds.
        cvNumFolds = 10;
        cvOptions = statset('crossval');
        regressMethod = @ReducedRankRegress;
        cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
            (regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
            numDimsUsedForPrediction, 'LossMeasure', 'NSE','RidgeInit', ...
            false, 'Scale', false);
        
        % when the partition is three-ways, j==1 means same target/source
        % pair and j==2 means diff target/source pair
        for j = [HPC, PFC]
            %             disp("processing rrr for "+p+" partition and the "+i+" pattern "+j+" direction")
            curr_source = (Patterns(p,j,i).X_source)';
            curr_target = (Patterns(p,j,i).X_target)';
            [   Patterns(p,j,i).rankRegress.cvl, ...
                Patterns(p,j,i).rankRegress.cvLoss,...
                Patterns(p,j,i).rankRegress.optDimReducedRankRegress,...
                Patterns(p,j,i).rankRegress.B,...
                Patterns(p,j,i).rankRegress.B_,...
                Patterns(p,j,i).rankRegress.V] ...
                = rankRegressRoutine(cvFun, cvNumFolds, ...
                cvOptions, curr_target, curr_source, ...
                numDimsUsedForPrediction);
            
            % Single neuron prediction
            for k = 1:nSource
                curr_singlesource = curr_source(:,j);
                if clean.zeroFiring(curr_singlesource)
                    continue;
                end
                [~,~, ...
                    dim_singleprediction{k}, ...
                    B_singleprediction{k},~,~] = ...
                    rankRegressRoutine(cvFun, cvNumFolds, ...
                    cvOptions,curr_target, ...
                    curr_singlesource,...
                    numDimsUsedForPrediction);
            end
            Patterns(p,j,i).rankRegress.singlesource_B = B_singleprediction;
            Patterns(p,j,i).rankRegress.singlesource_optDim = ...
                dim_singleprediction;
            Patterns(p,j,i).rankRegress.B_rrr = getReducedB_(Patterns(p,j,i).rankRegress.B,...
                Patterns(p,j,i).rankRegress.V, nSource, nTarget,...
                Patterns(p,j,i).rankRegress.optDimReducedRankRegress);
        end
    end
end


%%
%%%%%%%%%%%%%%%% FACTOR ANALYSIS SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doFactorAnalysis = false;
if doFactorAnalysis
    q = 1: nSource;
    cvNumFolds = 10;
    cvOptions = statset('crossval');
    regressMethod = @FactorRegress;
    for p = 1:Option.numPartition
        for i = 1:nPatternAndControl
            for j = [HPC, PFC]
                disp("processing the "+p+" partition and the "+i+" pattern "+j)
                currSource = Patterns(p,j,i).X_source';
                currTarget = Patterns(p,j,i).X_target';
                warning('')
                
                cvLoss = CrossValFa(currSource, q, cvNumFolds, cvOptions);
                qOpt =  FactorAnalysisModelSelect(cvLoss, q);
                
                cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
                    (regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
                    numDimsUsedForPrediction, ...
                    'LossMeasure', 'NSE', 'qOpt', qOpt);
                
                cvl = crossval(cvFun, ...
                    currTarget, currSource, ...
                    'KFold', cvNumFolds, ...
                    'Options', cvOptions);
                
                cvLoss = ...
                    [ mean(cvl); std(cvl)/sqrt(cvNumFolds) ];
                
                Patterns(p,j,i).factorAnalysis.optDimFactorRegress = ModelSelect...
                    (cvLoss, numDimsUsedForPrediction);
                
                Patterns(p,j,i).factorAnalysis.qOpt = qOpt;
                Patterns(p,j,i).factorAnalysis.cvl = cvl;
                Patterns(p,j,i).factorAnalysis.cvLoss = cvLoss;
                
                
                [Patterns(p,j,i).factorAnalysis.Z, ...
                    Patterns(p,j,i).factorAnalysis.U,...
                    Patterns(p,j,i).factorAnalysis.Q] = ...
                    ExtractFaLatents(currSource, qOpt);
                [warnMsg, warnId] = lastwarn();
                
                if isempty(warnId)
                    Patterns(p,j,i).singularWarning = false;
                    
                else
                    Patterns(p,j,i).singularWarning = true;
                    disp("singular is marked")
                end
            end
        end
    end
end

%%

for p = 1:Option.numPartition
    for i = 1:nPatternAndControl
        for j = [HPC, PFC]
            try
                [Patterns(p,j,i).rankRegress.removedPerformance, ~,~]=...
                    plots.getUncorrelatedPerformance( Patterns(p,j,i).rankRegress.B_,...
                    Patterns(p,j,i).X_source, Patterns(p,j,i).X_target, Patterns(p,j,i).rankRegress.optDimReducedRankRegress,...
                    numDimsUsedForPrediction, Patterns(p,j,i).rankRegress.cvLoss);
            catch
                keyboard
            end
            [Patterns(p,j,i).perfAsRemovingEach] = plots.sequentialRemovePredDims(Patterns(p,j,i).X_source, ...
                Patterns(p,j,i).X_target,Patterns(p,j,i).rankRegress.B_, Patterns(p,j,i).rankRegress.optDimReducedRankRegress, ...
                Patterns(p,j,i).rankRegress.cvLoss,numDimsUsedForPrediction);
        end
    end
end

%%

%%%%%%%%%%%%%%%% CREATE TABLE AND SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare table!
path2 = "C:\Users\BrainMaker\MATLAB Drive\Shared";
cd(path2)
if exist("Megatable.mat", 'file')
    load("Megatable.mat");
else
    TABLE = [];
end

if exist("T.mat", 'file')
    load("T.mat");
else
    warning("no exisiting real table")
    T = [];
end

params.optimizeOptionsScript;
Optiontable  = struct2table(Option, 'AsArray', true);
Patterntable = query.getPatternTable(Patterns);

if sum(Patterntable.singularWarning) > 0.7*size(Patterntable,1)
    warning("too many singular results in factor analysis")
    choice = input ("proceed to store results?");
    if choice == "no"
        error ("unsucessful run")
    end
end

% Identifying information about this options set and date of run
hash = DataHash(Option);
hash = hash(1:7); % Take the first 7 letters of the hash
hash = string(hash);
timestamp = string(date());
tablerow = [Optiontable, table(timestamp, hash, numWindowsCut, cutoffs, optimizationResult)]; % Combine option columnns with hash and date
tablecontent = [Patterntable, repmat(tablerow, height(Patterntable), 1)]; % combine those with all rows of the pattern table

%% Check and Hash

if ~isempty('TABLE') && any(contains(TABLE.hash, hash))
    TABLE(contains(TABLE.hash, hash), :) = []; % Delete any rows that contain the current hash
    T(contains(T.hash, hash), :) = [];
    disp("already computed before, rehashing to the same location");
    % New options:    Append row
else
    disp("new results stored!")
end

old_height = height(T);


% append the new results
if (size(T,2) ~= size(tablecontent,2))
    T = table.addNewColumn(T, tablecontent);
    TABLE    = table.addNewColumn(TABLE, tablerow);
else
    T = [T; tablecontent];
    TABLE = [TABLE; tablerow];
end
assert(height(T) > old_height, "not appending");

%% Save
% control.comparOptDim
cellInfo.getCellIdentities;
save ("Megatable", "TABLE",'-v7.3');
save ("T", "T", '-v7.3');
% writetable(TABLE,'~/Megatable.xls');
save(fullfile(datadefine, "hash", hash), ...
    "Patterns", "Option", "H", "cellOfWindows", "tablecontent","celllookup",'-v7.3')


% most recent state
save(fullfile(datadefine, 'mostRecentState'))
