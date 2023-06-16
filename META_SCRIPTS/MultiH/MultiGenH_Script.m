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
Default.animal    = "JS14";
%Default.generateH = cell(1,2);
Default.generateH{1} = "fromFilteredEEG "+" fromRipTimes";
Default.generateH{2} = "fromCoherence "+" fromRipTimes";
Default.generateH{3} = "fromWpli "+" fromRipTimes";
Default.samplingRate  = [] ;              % For spikes.getSpikeTrain, nan if not given
Default.spikeBinSize  = 0.15;               % 100 milliseconds
Default.timesPerTrial = 10;                % 10 times per trial
Default.winSize       = {[-0.15, 0.15]};     % size of the window
Default.sourceArea    = "CA1";             % only when there are
Default.equalWindowsAcrossPatterns = true;    % whether all three patterns have the same #windows
Default.singleControl = false;                 % whether to use just one control column
Default.numPartition = 50;                    % ways to split source and target
Default.binsToMatchFR = 20;
%Default.binsToMatchWindows = 50;
Default.quantileToMakeWindows = 0.85;
Default.preProcess_FilterLowFR = true;
Default.preProcess_matchingDiscreteFR = true;
Default.shortedRippleWindow = false; % using 100 milisecond window for ripples
Default.oldControlBehavior = false;
Default.lowerControl = true;

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
cellOfWindows = cell(numel(Option.generateH),1);
frequenciesPerPattern = [6 10; 0.5 4;  150 200];
[nPatterns,~] = size(frequenciesPerPattern);

if Option.singleControl == true
    nPatternAndControl = nPatterns+1;
else
    nPatternAndControl = nPatterns*2;
end

X_pfc = cell(numel(Option.generateH),1);
X_hpc = cell(numel(Option.generateH),1);

%%%%%%%%%%%%%%%% READ IN RAW PATTERN DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cutoffs = zeros(numel(Option.generateH), nPatterns*2);
Hc_cutoffs = zeros(numel(Option.generateH), nPatterns*2);

numWindowsCut = zeros(numel(Option.generateH), nPatterns);

%% LOOP THROUGH THE GEN-H METHODS
% -------------------------------
for i = 1:numel(Option.generateH)
    
    currMethod = Option.generateH{i};
    
    % Find network pattern events
    if contains(currMethod, "fromSpectra")
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
    elseif contains(currMethod, "fromFilteredEEG")
        load(Option.animal + "avgeeg.mat");
         [~, sleepSessions] = getRunningSessions(Option.animal);
        [H, Hvals, Hnanlocs, Htimes] = eventMatrix.generateFromFilteredEEG(avgeeg, ...
            Option.sourceArea, "patterns",patternNames(1:3),"downsample",10, "sleepSessions", sleepSessions);
    elseif contains(currMethod, "fromCoherence")
        load(Option.animal + "spectralBehavior.mat");
        spectrogram = efizz.C;
        frequencyAxis = efizz.f;
        Htimes = efizz.t;
        [H, Hvals, Hnanlocs, Htimes] = eventMatrix.generateFromSpectra(Htimes, spectrogram, frequencyAxis,...
            frequenciesPerPattern);
    elseif contains(currMethod, "fromWpli")
        load(Option.animal + "spectralBehavior.mat");
        spectrogram = efizz.wpli;
        frequencyAxis = efizz.f;
        Htimes = efizz.t;
        [H, Hvals, Hnanlocs, Htimes] = eventMatrix.generateFromSpectra(Htimes, spectrogram, frequencyAxis,...
            frequenciesPerPattern);
        
    else
        error("Core method for deriving the event matrix is not recognized");
    end
    
    %% getting correlations
    corrMatrix = params.getPatternSpeedCorr(animal, H, Hvals, Htimes);
    theta_sp(i) = corrMatrix(2,1);
    theta_ripple(i) = corrMatrix(4,2);
    
    %% Modify ripple pattern? Lower the threshold as with window sizes
    
    if contains(currMethod,"fromRipTimes")
        
        load(Option.animal + "globalripple01.mat");
        % COHERENCE/WPLI RIPPLE
        if any(contains(Option.generateH{i}, ["fromWpli", "fromCoherence"]))
            [~, H(:,RIPPLE), Hnanlocs(:,RIPPLE), Hvals(:,RIPPLE), minRippleThreshold, original] = ...
                eventMatrix.generateFromRipples(globalripple, ...
                'amplitude_at_riptime', true,...
                'rippleBand', Hvals(:,RIPPLE),...                               RY: Hvals, not H here for obvious reasons: you  want the original ripple band activity
                'rippleBandTime', Htimes,...
                'globalrippleWindowUnits', 'amp');
            % HILBERT/FFT RIPPLE
        else
            [~, H(:,RIPPLE), Hnanlocs(:,RIPPLE), Hvals(:,RIPPLE), minRippleThreshold, original] = ...
                eventMatrix.generateFromRipples(globalripple, ...
                'amplitude_at_riptime', true,...
                'rippleBand', Hvals(:,RIPPLE),...                               RY: Hvals, not H here for obvious reasons: you  want the original ripple band activity
                'rippleBandTime', Htimes,...
                'globalrippleWindowUnits', 'std');
        end
        
    end
    
    %% Making Windows
    %%%%%%%%%%%%%%%% WINDOW SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [cellOfWindows{i}, cutoffs(i,THETA:DELTA)] = windows.make(   Htimes, Option.quantileToMakeWindows, H(:,THETA:DELTA), winSize);
    if any(contains(Option.generateH{i}, ["fromCoherence","fromWpli"]))
        [cellOfWindows{i}(RIPPLE), cutoffs(i,RIPPLE)] = windows.make(Htimes,...
            Option.quantileToMakeWindows,  H(:,RIPPLE), winSize,...
            'quantile', Hvals(:,RIPPLE),'higherThanQuantile', true); % % RY: quantile needs to be hvals for ripple coherence/wpli threshold to be correct, but timesd computed from H such that non-ripple times thrown out
    else
        [cellOfWindows{i}(RIPPLE), cutoffs(i,RIPPLE)] = windows.make(Htimes, ...
            1,   H(:,RIPPLE), winSize, 'threshold', 'raw','higherThanQuantile', true);
    end
    disp(newline);
    disp('Initiale window  creation:')
    disp('--------------------------')
    cellfun(@(x,y)  fprintf("%d windows for %s\n", size(x,1),  y), ...
        cellOfWindows{i}, cellstr(patternNames(1:3)));
    
    % equalize number of windows across patterns based on input argument
    if Option.equalWindowsAcrossPatterns == true
        cellOfWindows{i} = windows.equalizeWindowsAcrossPatterns(cellOfWindows{i}, nPatterns);
    end
    
    disp(newline);
    disp('Equalized windows        :')
    disp('--------------------------')
    cellfun(@(x,y)  fprintf("%d windows for %s\n", size(x,1),  y), ...
        cellOfWindows{i}, cellstr(patternNames(1:3)));
    
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
    
    % Old way, scramble and randomly select times
    if Option.oldControlBehavior
        Hc =  control.generatePatternShuffle(H(:,1:3), Htimes, cellOfWindows{i}); % add control patterns;
        % New way, low activity!
    else
        Hc = H;
    end
    
    % Theta/Delta
    [Hc_cellOfWindows{i}, Hc_cutoffs(i,THETA:DELTA)] = windows.make(Htimes,  quantileControl,...  % add windows of control patterns
        Hc(:,THETA:DELTA), winSize,... % Selects less than quantile
        'higherThanQuantile', Option.oldControlBehavior);
    
    % Ripples
    % -------
    % LOW RIPPLE WPLI/COHERENCE
    if any(contains(Option.generateH{i}, ["fromCoherence","fromWpli"]))
        [Hc_cellOfWindows{i}(RIPPLE), Hc_cutoffs(i,RIPPLE)] = windows.make(Htimes,...
            quantileControl,  Hvals(:,RIPPLE), winSize,... 01/03/2020 Hc should be Hvals so as to not have nans in non-ripple periods
            'quantile', Hvals(:,RIPPLE),'higherThanQuantile', Option.oldControlBehavior);
        % LOW RIPPLE FFT/HILBERT
    else
        if Option.oldControlBehavior
            [Hc_cellOfWindows{i}(RIPPLE), Hc_cutoffs(i,RIPPLE)] = windows.make(Htimes, ...
                1,   Hc(:,RIPPLE), winSize, 'threshold', 'raw','higherThanQuantile', true); %RY quantile won't work because these are raw
        else
            [Hc_cellOfWindows{i}(RIPPLE), Hc_cutoffs(i,RIPPLE)] = windows.make(Htimes, ...
                quantileControl,   Hvals(:,RIPPLE), winSize, 'threshold', 'quantile','higherThanQuantile', false); %RY quantile won't work because these are raw
        end
    end
    
    % Clean
    % -----
    % clean up control windows: remove each control pattern's window's overlap
    for pattern = 1:nPatterns
        curr = windows.removeOverlapsBetweenPattern(...
            cell2mat(cellOfWindows{i}(:,pattern)), cell2mat(Hc_cellOfWindows{i}(:,pattern)));
        Hc_cellOfWindows{i}{pattern} = curr;
    end
    
    % Merge many controls into 1 control?
    % -----------------------------------
    % % Merge into one
    cellOfWindows{i}(nPatterns+1:nPatterns*2) = Hc_cellOfWindows{i};
    %     cutoffs = [cutoffs,Hc_cutoffs];
    
    % Ensure each pattern has equal # of window
    % -----------------------------------------
    % Equalize trials/windows for each pair of patttern-controlPattern
    [cellOfWindows{i}, warnedEmptyControls] =...
        control.equalizePatternControl(nPatterns, cellOfWindows{i});
    
    
    
    if ~any(contains(patternNames,"control"))
        patternNames = [patternNames; patternNames+"-control"]';
        patternNames = patternNames(:)';
    end
    
    
    %% Getting spikes
    [timeBinStartEnd, timeBinMidPoints, times_spiking, spikeCountMatrix, spikeRateMatrix, ...
    areaPerNeuron, cell_index, sessionTypePerBin] = spikes.getSpikeTrain(Option.animal, ...
        Option.spikeBinSize,  Option.samplingRate);
    
    if Option.preProcess_FilterLowFR % filter the neurons whose firing rate is lower than specified threshold
        [spikeCountMatrix, spikeRateMatrix, avgFR, areaPerNeuron, cell_index]...
            = trialSpikes.filterFR(spikeCountMatrix, spikeRateMatrix, 0.1, ...
            timeBinStartEnd, areaPerNeuron, cell_index);
    end
    [spikeSampleMatrix, spikeSampleTensor, trialTimes] = trialSpikes.generate(spikeCountMatrix,...
        timeBinStartEnd, cellOfWindows{i}, ...
        Option.timesPerTrial, nPatternAndControl);
     arrayfun(@(i) range(trialTimes{i}(:))/3600,1:6)
    %% Separate spikesSampleMatrix/Tensor by area that neurons are in PFC and neurons that in HPC
   
    [pfcFR, hpcFR] = trialSpikes.separateFiringRate(avgFR, areaPerNeuron);
    X_pfc{i} = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "PFC");
    X_hpc{i} = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "CA1");
    
    %% Separate firing pattern into source and target
    [nPFCneurons,~] = size(X_pfc{i}{1});
    [nHPCneurons,~] = size(X_hpc{i}{1});
end
%%
%%%%%%%%%%%% Assign all of the RAW relevent structures %%%%%%%%%%%%
Raw = struct();
kws = {'UniformOutput',false};
Raw.X_pfc = cellfun(@(x) cellfun(@(y) single(y), x, kws{:}), X_pfc, kws{:});
Raw.X_hpc = cellfun(@(x) cellfun(@(y) single(y), x, kws{:}), X_hpc, kws{:});
Raw.H = struct('H', H,...
    'Hvals', Hvals,...
    'Hnanlocs', Hnanlocs,...
    'Htimes', Htimes)
Raw.frequenciesPerPattern = frequenciesPerPattern;
Raw.cellOfWindows = cellOfWindows;
Raw.spikeSampleTensor = spikeSampleTensor;

% -------------------- Generate empty pattern struct --------------------
%% Results place to store outputs
clear Patterns
% Patterns become a cell of two struct arrays?
% Patterns = struct();

Patterns = struct("X_source",[], "X_target",[]);
Patterns.rankRegress = struct(...
    "B", [], ...
    "B_", [], ...
    "optDimReducedRankRegress", 0, ...
    "singlesource_B", [], ...
    "singlesource_optDim",[]);
Patterns.factorAnalysis = struct("qOpt", []);

patternNames = ["theta","delta","ripple",...
    "theta-control","delta-control","ripple-control"];

if nPatternAndControl == nPatterns+1
    patternNames = ["theta","delta","ripple","control"];
end

Patterns = repmat(Patterns, ...
    [numel(Option.generateH), Option.numPartition, 2, nPatternAndControl]);

if Option.sourceArea == "CA1"
    directionality = ["hpc-hpc","hpc-pfc"];
else
    directionality = ["pfc-hpc","pfc-pfc"];
end

% -------------------- Assign cell partitions per pattern ---------------

% Label all of the non-cell related dimensions of the pattern
Patterns = nd.dimLabel(Patterns, 1, "generateH",      string(Option.generateH)); % Label the generateH dimension of struct
Patterns = nd.dimLabel(Patterns, 2, "iPartition"); % Label the generateH dimension of struct
Patterns = nd.dimLabel(Patterns, 3, "directionality", directionality);           % Label the directionality dimension of struct
Patterns = nd.dimLabel(Patterns, 4, "name",           patternNames);             % Label the directionality dimension of struct
%%
% Calculate and place the cell-wise partition splits
for iPartition = 1:Option.numPartition
    
    if lower(Option.sourceArea) == "ca1"
        s_hpc = cell(numel(Option.generateH),1); % source hpc for two methods
        t_hpc = cell(numel(Option.generateH),1); % target hpc for two methods
        t_pfc = cell(numel(Option.generateH),1); % target pfc for two methods
        
        % the first method
        [s_hpc{1}, target, nSource, nTarget,...
            s_hpc_index, index_target] = ...
            trialSpikes.matchFRinDiscreteRanges...
            (X_hpc{1}, X_pfc{1}, hpcFR, pfcFR, Option.binsToMatchFR);
        
        t_hpc{1} = target(1,:);
        t_pfc{1} = target(2,:);
        t_hpc_index = index_target(1,:);
        t_pfc_index = index_target(2,:);
        
        % the second method
        for j = 2:numel(Option.generateH)
            for i = 1:numel(patternNames)
                s_hpc{j}{i} = X_hpc{j}{i}(s_hpc_index,:);
                t_hpc{j}{i} = X_hpc{j}{i}(t_hpc_index,:);
                t_pfc{j}{i} = X_pfc{j}{i}(t_pfc_index,:);
            end
        end
        
        % Ensure different
        for i = numel(patternNames)
            for j = numel(t_hpc)
                assert( ~isequal( t_hpc{j}{i}, t_pfc{j}{i} ) );
            end
        end
        
        %     elseif lower(Option.sourceArea) == "pfc"
        %         [s_pfc{1},target{1}, nSource, nTarget,...
        %             s_pfc, index_target] = ...
        %             trialSpikes.matchFRinDiscreteRanges...
        %             (X_hpc{1}, X_pfc{1}, hpcFR, pfcFR, Option.binsToMatchFR);
        %         [s_pfc{2},target{2}, nSource, nTarget,...
        %             s_pfc, index_target] = ...
        %             trialSpikes.matchFRinDiscreteRanges...
        %             (X_hpc{2}, X_pfc{2}, hpcFR, pfcFR, Option.binsToMatchFR);
        %         t_pfc = target{1}(1,:);
        %         t_hpc = target{1}(2,:);
        %         t_pfc = index_target(1,:);
        %         t_hpc = index_target(2,:);
    end
    % ------------------------------
    % Place paritioned data properly
    % ------------------------------
    for k = 1:numel(Option.generateH)
        for i = 1:numel(patternNames)
            for j = 1:numel(directionality)
                % Parse directionality
                
                sourcetarg = directionality(j).split('-');
                
                source = sourcetarg(1); % current source region
                target = sourcetarg(2); % current target region
                
                if source == "hpc"
                    s_dat = s_hpc{k}{i};
                    s_ind = s_hpc_index;
                elseif source == "pfc"
                    s_dat = s_pfc{k}{i};
                    s_ind = s_pfc_index;
                end
                if target == "hpc"
                    t_dat = t_hpc{k}{i};
                    t_ind = t_hpc_index;
                elseif target == "pfc"
                    t_dat = t_pfc{k}{i};
                    t_ind = t_pfc_index;
                end
                
                % Assign x_source and x_target
                Patterns(k,iPartition,j,i).X_source = single(s_dat);
                Patterns(k,iPartition,j,i).X_target = single(t_dat);
                
                % Assign index_source and index_target
                Patterns(k,iPartition,j,i).index_source = s_ind;
                Patterns(k,iPartition,j,i).index_target = t_ind;
            end
        end
    end
end

% Ensure different
indices = nd.indicesMatrixForm(Patterns);
indices = unique(indices(:,[1 2 4]), 'rows');
for index = indices'
    I = num2cell(index);
    A = Patterns(I{1:2}, 1, I{3}).X_target;
    B = Patterns(I{1:2}, 2, I{3}).X_target;
    assert( ~isequal(A, B) );
end

Raw.patternNames = patternNames;
Raw.directionality = directionality;

%%
%%%%%%%%%%%%%%%% RANK-REGRESS SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nTarget = size(Patterns(1,1,1,1).X_target,1);
nSource = min(size(Patterns(1,1,1,1).X_source,1),...
    size(Patterns(1,1,2,1).X_source,1));

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
        for k = 1:numel(Option.generateH)
            % when the partition is three-ways, j==1 means same target/source
            % pair and j==2 means diff target/source pair
            for j = [HPC, PFC]
                
                %             disp("processing rrr for "+p+" partition and the "+i+" pattern "+j+" direction")
                curr_source = double((Patterns(k,p,j,i).X_source)');
                curr_target = double((Patterns(k,p,j,i).X_target)');
                [   Patterns(k,p,j,i).rankRegress.cvl, ...
                    Patterns(k,p,j,i).rankRegress.cvLoss,...
                    Patterns(k,p,j,i).rankRegress.optDimReducedRankRegress,...
                    Patterns(k,p,j,i).rankRegress.B,...
                    Patterns(k,p,j,i).rankRegress.B_,...
                    Patterns(k,p,j,i).rankRegress.V] ...
                    = rankRegressRoutine(cvFun, cvNumFolds, ...
                    cvOptions, curr_target, curr_source, ...
                    numDimsUsedForPrediction);
                
                % Single neuron prediction
                for h = 1:nSource
                    curr_singlesource = curr_source(:,j);
                    if clean.zeroFiring(curr_singlesource)
                        continue;
                    end
                    [~,~, ...
                        dim_singleprediction{h}, ...
                        B_singleprediction{h},~,~] = ...
                        rankRegressRoutine(cvFun, cvNumFolds, ...
                        cvOptions,curr_target, ...
                        curr_singlesource,...
                        numDimsUsedForPrediction);
                end
                Patterns(k,p,j,i).rankRegress.singlesource_B = B_singleprediction;
                Patterns(k,p,j,i).rankRegress.singlesource_optDim = ...
                    dim_singleprediction;
                Patterns(k,p,j,i).rankRegress.B_rrr = getReducedB_(Patterns(k,p,j,i).rankRegress.B,...
                    Patterns(k,p,j,i).rankRegress.V, nSource, nTarget,...
                    Patterns(k,p,j,i).rankRegress.optDimReducedRankRegress);
            end
        end
    end
end


%%%%%%%%%%%%%%%% SAVE RESULTS OF PATTERNS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --------------
%% PART 1 : TABLE
%% --------------
 Patterntable = [];
optionRun = struct2table(Option, 'AsArray', true);

% unique for each option run
hash = DataHash(Option);
hash = hash(1:8);
hash = string(hash);
timestamp = string(date());


for k = 1:numel(Option.generateH)
    OptionToStore = rmfield(Option, 'generateH');
    OptionToStore.generateH = Option.generateH{k};
    Optiontable = struct2table(OptionToStore, 'AsArray', true);
    
    [optResult, perf] = params.optmizeOptions...
        (theta_sp(k), theta_ripple(k), Patterns(k,:,:,:), nTarget, nSource, Option.numPartition);
    currPatterntable  = query.getPatternTable(Patterns(k,:,:,:),Option);
    windowCutOff      = cutoffs(k,:);
    theta_sp_corr     = theta_sp(k);
    theta_ripple_corr = theta_ripple(k);
    optResultsrow = table(windowCutOff,theta_sp_corr, theta_ripple_corr ,...
        perf, optResult);
    
    tablerow = [table(timestamp, hash), optResultsrow]; % Combine option columnns with hash and date
    
    currPatterntable.generateH = [];
    currPatterntable = [currPatterntable, repmat(tablerow,height(currPatterntable),1)];
    Patterntable = [Patterntable; currPatterntable];
    
end

%% -------------------------
%% PART 2 : SAVE HASHED DATA
%% -------------------------
% most recent state
% save(fullfile(datadefine, 'mostRecentState'))

%% -------------------
%% PART 3 : SAVE TABLE
%% -------------------
path2 = "C:\Users\BrainMaker\MATLAB Drive\Shared";
cd(path2)
% if exist("TABLE_multi3.mat", 'file')
%     load("TABLE_multi3.mat");
%     load("Tablerow_3.mat");
% else
%     TABLE_multi = [];
%     Tablerow = [];
% end

if exist("TABLE_multi3_FRfixed.mat", 'file')
    load("TABLE_multi3_FRfixed.mat");
    load("Tablerow_3_FRfixed.mat");
else
    TABLE_multi3_FRfixed = [];
    Tablerow_3_FRfixed = [];
end
%optionRun = [optionRun, table(convertCharsToStrings(hash))];
%optionRun.Properties.VariableNames{18} = 'hash'; % No bueno : if optionRun changes column size, this will break. Worse, if I cut-paste/swapping lines that create the option struct at the beginning, this will also break it
optionRun.hash = string(hash);
% append the new results

if size(TABLE_multi3_FRfixed,2) ~= size(Patterntable,2)
    Tablerow = table.addNewColumn(Tablerow, optionRun);
    TABLE_multi    = table.addNewColumn(TABLE_multi, Patterntable);
else
    Tablerow_3_FRfixed = [Tablerow_3_FRfixed; optionRun];
    TABLE_multi3_FRfixed = [TABLE_multi3_FRfixed; Patterntable];
end

%%
save ("TABLE_multi3_FRfixed", "TABLE_multi3_FRfixed",'-v7.3');
save ("Tablerow_3_FRfixed", "Tablerow_3_FRfixed", '-v7.3');

save(fullfile(datadefine, "hash_FRfixed", hash), ...
    "Patterns", "Option", "Raw", '-v7.3') % 3 structs that carry everything needed about the run
