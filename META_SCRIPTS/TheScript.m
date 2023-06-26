% TODO:
% 1. have way of passing generate_H cell, and then just load all of those
%    patterns simulateously in 1 go instead of as separate runs
% 
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
disp("------------------------")
disp("    Obtaining events    ")
disp("------------------------")
[Events] = events.ThetaDeltaRipple(Option);
% Documentation
% Events is a struct with fields:
% - .times : array of times of events
% - .H     : Event Matrix,    T x 3, and each column are theta, delta, ripple
% - .Hvals : Event Matrix,    T x 3, values without nans
% - .Hnanlocs : Event Matrix, T x 3, logicals of nans

%%%%%%%%%%%%%%%% CUT WINDOWS WITH EVENT MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("------------------------")
disp("    Cutting windows     ")
disp("------------------------")
[cellOfWindows, cutoffs] = windows.ThetaDeltaRipple(Events, Option);
% Documentation
% -  cellOfWindows: 1 x nPatterns cell array of windows
% -  cutoffs:       nPatterns x 1 vector of cutoffs
% TODO: modify to be able to include overall pattern and track patterns
% PRIORITY; overall: medium, track: very low, overall can be included in
% cellOfWindows, whereas, track can be included as a separate output

%%%%%%%%%%%%%%%% ACQUIRE SPIKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Getting spikes
disp("------------------------")
disp("    Getting spikes      ")
disp("------------------------")
[timeBinStartEnd, timeBinMidPoints, ~, spikeCountMatrix, spikeRateMatrix, ...
    areaPerNeuron, cell_index, sessionTypePerBin] = ...
    spikes.getSpikeTrain(Option.animal, ...
            Option.spikeBinSize,  Option.samplingRate);

% filter the neurons whose firing rate is lower than specified threshold
if Option.preProcess_FilterLowFR 
    disp("------------------------")
    disp("Filtering low FR neurons")
    disp("------------------------")
    [spikeCountMatrix, spikeRateMatrix, avgFR, areaPerNeuron, cell_index]...
        = trialSpikes.filterFR(spikeCountMatrix, spikeRateMatrix, 0.1, ...
        timeBinStartEnd, areaPerNeuron, cell_index);
end


%%%%%%%%%%%%%%%% ACQUIRE TRIALS FROM WINDOWS + SPIKES %%%%%%%%%%%%%%%%%%%
% RYAN bug here .. timeBinStartEnd instead of timeBinMidPoints
disp("------------------------")
disp("   Windowing spikes     ")
disp("------------------------")
[spikeSampleMatrix, spikeSampleTensor, trialTimes] = trialSpikes.generate(...
    spikeCountMatrix, timeBinMidPoints, cellOfWindows, ... 
    Option.timesPerTrial, Option.nPatternAndControl);


% %%%%%%%%%%%%%%% SETUP RAW DATA STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure for separated data
r.hpc = struct;
r.pfc = struct;
%%%%%%%%%%%%%%%% SEPRATE BRAIN AREA DATASULT STRUCTURES %%%%%%%%%%%%%%%%%%
% Separate spikesSampleMatrix/Tensor by area that neurons are in PFC and
% neurons that in HPC
[r.pfc.FR, r.hpc.FR] = trialSpikes.separateFiringRate(avgFR, areaPerNeuron);
r.pfc.T = trialSpikes.separateSpikes(spikeSampleTensor, areaPerNeuron, "PFC");
r.hpc.T = trialSpikes.separateSpikes(spikeSampleTensor, areaPerNeuron, "CA1");
r.pfc.X = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "PFC");
r.hpc.X = trialSpikes.separateSpikes(spikeSampleMatrix, areaPerNeuron, "CA1");
r.trialTimes = trialTimes;

%% Separate firing pattern into source and target
[nPFCneurons,~,~] = size(r.pfc.X{1});
[nHPCneurons,~,~] = size(r.hpc.X{1});
r.celllookup = cellInfo.getCellIdentities(Option.animal, cell_index,...
                                          areaPerNeuron);
r.avgFR         = avgFR;
r.areaPerNeuron = areaPerNeuron;
r.pfc.nNeurons  = nPFCneurons;
r.hpc.nNeurons  = nHPCneurons;
r.windowInfo.cellOfWindows = cellOfWindows;
r.windowInfo.nWindows      = cellfun(@(x) size(x, 1), cellOfWindows);
r.nPattern = Option.nPatterns;
r.nControl = Option.nPatternAndControl - r.nPattern;
r.timeBinMidPoints = timeBinMidPoints;
r.sessionTypePerBin = sessionTypePerBin;
r.spikeRateMatrix = spikeRateMatrix;
r.spikeCountMatrix = spikeCountMatrix;

%%%%%%%%%%%%%%%% SETUP PARTITIONS AND RESULT STRUCTURES %%%%%%%%%%%%%%%%%%
disp("------------------------")
disp(" Subsampling partitions ")
disp("------------------------")
[Patterns, Patterns_overall] = trialSpikes.partitionAndInitialize(r, Option);

%%%%%%%%%%%%%%%% ANALYSIS SECTION    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("------------------------")
disp("     Analysis           ")
disp("------------------------")


if Option.analysis.rankRegress
    % Rank regression of network pattern windows of spiking activity
    % (Subspaces acquired here)
    % TODO: 
    % 1. fix Option.rankregress => Option.rankRegress                 
    % 2. most rankRegress.B_ are empty                                
    Patterns         = analysis.rankRegress(Patterns, Option);        
    Patterns_overall = analysis.rankRegress(Patterns_overall, Option);
end

if Option.analysis.factorAnalysis
    % Factor analysis of network pattern windows of spiking activity
    % (Used to measure instrinsic dimensionality of network activity)
    Patterns = analysis.factorAnalysis(Patterns, Option);
end

if Option.analysis.cca
    % BUG: mismatch in size of source matrix and CCA u v matrices
    Patterns         = analysis.cca(Patterns, Option);
    Patterns_overall = analysis.cca(Patterns_overall, Option);
    % TODO : section that knocks off kim 2022 after these measurements
end

if Option.analysis.timeVarying
    % How much spiking moment to moment is explained by subspace
    % ISSUE: hits a bug on line 4
    % TODO: 1 .also return epochwise zscored neural firing matching
    %       2. return timeseries of smoothed firing rate
    running_times = r.timeBinMidPoints(r.sessionTypePerBin == 1);
    [behavior, thrown_out_times] = table.behavior.lookup(Option.animal, ...
                                                         running_times);
    Components = analysis.timeVarying_v2(Patterns, Option, r);
    Components = plots.temporal.correlateSpectral(Components, Events, Option);
    Components = plots.temporal.correlateBehavior(Components, Events, Option);
end

if Option.analysis.checks
    % Plots regarding the raw and processed data (and sometimes
    % relation to processed Patterns struct)
    % TODO: Think about splitting this into checks involving
    %        versus not involving the Patterns struct
    plots.runChecks(Patterns, Option, r, cellOfWindows);
end

% TODO: (1) plug in JPECC version of rankRegress here
% TODO: (2) function that outputs average response of Pattern struct per neuron

%%%%%%%%%%%%%%%% CREATE TABLE AND SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Option.save

    disp("Saving results")

    % Ready the path
    table_folder = fullfile(codedefine(), 'DATA_TABLES');
    if ~exist(table_folder, 'dir')
        mkdir(table_folder);
    end
    if ~ismember(table_folder, path)
        addpath(table_folder);
    end

    %% RunsSummary: Abbreviated summary of results for runs
    if exist("RunsSummary.mat", 'file') 
        load("RunsSummary.mat");
    else
        RunsSummary = [];
    end

    %% DetailedRunsSummary: Summary of information for runs
    if exist("DetailedRunsSummary.mat", 'file') 
        load("DetailedRunsSummary.mat");
    else
        warning("no existing real table")
        DetailedRunsSummary = [];
    end

    % Determine information to add to table
    params.optimizeOptionsScript;
    %Option = rmfield(Option, "animal");
    Optiontable  = struct2table(Option, 'AsArray', true);
    Patterntable = query.getPatternTable(Patterns);

    % Check for singular results in factor analysis
    factorAnalysisChecksum = sum(Patterntable.singularWarning) ...
                             > 0.7 * size(Patterntable,1);
    if factorAnalysisChecksum
        warning("too many singular results in factor analysis")
        choice = input ("proceed to store results?");
        if choice == "no"
            error ("unsuccessful run")
        end
    end

    % Identifying information about this options set and date of run
    hash = DataHash(Option);
    hash = hash(1:7); % Take the first 7 letters of the hash
    hash = string(hash);
    timestamp = string(date());
    tablerow = [Optiontable, ...
        table(timestamp, hash, numWindowsCut, cutoffs, optimizationResult)]; % Combine option columns with hash and date
    tablecontent = [Patterntable, repmat(tablerow, height(Patterntable), 1)]; % combine those with all rows of the pattern table

    %% Check and Hash
    if ~isempty('RunsSummary') && any(contains(RunsSummary.hash, hash))
        RunsSummary(contains(RunsSummary.hash, hash), :) = []; % Delete any rows that contain the current hash
        DetailedRunsSummary(contains(DetailedRunsSummary.hash, hash), :) = [];
        disp("already computed before, rehashing to the same location");
        % New options:    Append row
    else
        disp("new results stored!")
    end

    % append the new results
    old_height = height(DetailedRunsSummary);
    if size(DetailedRunsSummary,2) ~= size(tablecontent,2)
        DetailedRunsSummary = table.addNewColumn(DetailedRunsSummary, tablecontent);
        RunsSummary    = table.addNewColumn(RunsSummary, tablerow);
    else
        DetailedRunsSummary = [DetailedRunsSummary; tablecontent];
        RunsSummary = [RunsSummary; tablerow];
    end
    assert(height(DetailedRunsSummary) > old_height, "not appending");

    %% ------------- Save ----------------------------
    % save the tables
    save("RunsSummary", "RunsSummary",'-v7.3');
    save("DetailedRunsSummary", "DetailedRunsSummary", '-v7.3');
    % save the results
    save(fullfile(datadefine, "hash_forTesting", hash), ...
        "Option", "Behaviors","celllookup",'-v7.3')
    % most recent state
    save(fullfile(datadefine, 'mostRecentState'))
end

% TODO: consider whether this is needed ... should be r isntead?
if Option.saveRaw
    Raw = struct();
    kws = {'UniformOutput',false};
    Raw.X_pfc = cellfun(@(x) cellfun(@(y) single(y), x, kws{:}), X_pfc, kws{:});
    Raw.X_hpc = cellfun(@(x) cellfun(@(y) single(y), x, kws{:}), X_hpc, kws{:});
    Raw.H = struct('H', H,...
                   'Hvals', Hvals,...
                   'Hnanlocs', Hnanlocs,...
                   'Htimes', Htimes);
    Raw.frequenciesPerPattern = frequenciesPerPattern;
    Raw.cellOfWindows = cellOfWindows;
    Raw.spikeSampleTensor = spikeSampleTensor;
    Raw.patternNames = patternNames;
    Raw.directionality = directionality;
    Raw.celllookup = celllookup;
    save(fullfile(datadefine, "hash_Raw", hash), ...
        "Option", "Raw",'-v7.3')
end
