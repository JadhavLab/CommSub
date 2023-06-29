% TODO:
% 1. have way of passing generate_H cell, and then just load all of those
%    patterns simulateously in 1 go instead of as separate runs
% 
% ---- PATH -----
% Matlab uses startup.m to run startup code...
% Put this in your startup.m so that the code for this is in path:
%
% addpath(genpath('/Volumes/MATLAB-Drive/')) % or wherever your CODE files are
% located
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
Events = events.ThetaDeltaRipple(Option);
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
Events = windows.ThetaDeltaRipple(Events, Option);
% -  cutoffs:       nPatterns x 1 vector of cutoffs
% TODO: modify to be able to include overall pattern and track patterns
% PRIORITY; overall: medium, track: very low, overall can be included in
% cellOfWindows, whereas, track can be included as a separate output

%%%%%%%%%%%%%%%% ACQUIRE SPIKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Getting spikes
disp("------------------------")
disp("    Getting spikes      ")
disp("------------------------")
Spk = spikes.getSpikeTrain(Option.animal, Option.spikeBinSize, ...
                         Option.samplingRate);

% filter the neurons whose firing rate is lower than specified threshold
if Option.preProcess_FilterLowFR 
    disp("------------------------")
    disp("Filtering low FR neurons")
    disp("------------------------")
    Spk = trialSpikes.filterFR(Spk, 0.1);
    disp("Mean FR: " + sort(Spk.avgFR))
end

if Option.preProcess_gaussianFilter
    % Gaussian filter the spikeCountMatrix/spikeRateMatrix
    gauss = gausswin(Option.preProcess_gaussianFilter);
    for i = progress(1:size(Spk.spikeRateMatrix, 1), 'Title', 'Gaussian filtering')
        Spk.spikeRateMatrix(i, :)  = conv(Spk.spikeRateMatrix(i, :), gauss, 'same');
        Spk.spikeCountMatrix(i, :) = conv(Spk.spikeCountMatrix(i, :), gauss, 'same');
    end
end

if Option.preProcess_zscore
    % Z-score the spikeCountMatrix/spikeRateMatrix
    disp(" Z-scoring ")
    Spk.spikeRateMatrix  = zscore(Spk.spikeRateMatrix,  0, 2);
    Spk.spikeCountMatrix = zscore(Spk.spikeCountMatrix, 0, 2);
    Spk.avgFR = mean(Spk.spikeRateMatrix, 2);
end
prewindow_copy = Spk;


%%%%%%%%%%%%%%%% ACQUIRE TRIALS FROM WINDOWS + SPIKES %%%%%%%%%%%%%%%%%%%
% RYAN bug here .. timeBinStartEnd instead of timeBinMidPoints
disp("------------------------")
disp("   Windowing spikes     ")
disp("------------------------")
Spk = trialSpikes.generate(Spk, Events, Option);

% %%%%%%%%%%%%%%% SETUP RAW DATA STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure for separated data
%%%%%%%%%%%%%%%% SEPRATE BRAIN AREA DATASULT STRUCTURES %%%%%%%%%%%%%%%%%%
% Separate spikesSampleMatrix/Tensor by area that neurons are in PFC and
% neurons that in HPC
%% Separate firing pattern into source and target
[Spk.nSource,~,~] = size(Spk.hpc.X{1});
[Spk.nTarget,~,~] = size(Spk.pfc.X{1});
Spk.celllookup = cellInfo.getCellIdentities(Option.animal, Spk.cell_index,...
                                            Spk.areaPerNeuron);

%%%%%%%%%%%%%%%% SETUP PARTITIONS AND RESULT STRUCTURES %%%%%%%%%%%%%%%%%%
disp("------------------------")
disp(" Subsampling partitions ")
disp("------------------------")
[Patterns, Patterns_overall] = trialSpikes.partitionAndInitialize(Spk, Option);

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
    % ISSUE: warnings emitted regarding full rank -- linear independence
    % violation can be subtle problem or not a problem at all
    % remedies (1) regularize (2) remove linearly dependent columns (PCA)
    Patterns         = analysis.cca(Patterns, Option);
    Patterns_overall = analysis.cca(Patterns_overall, Option);
    % TODO : section that knocks off kim 2022 after these measurements
end

if Option.analysis.timeVarying
    % How much spiking moment to moment is explained by subspace
    % ISSUE: hits a bug on line 4
    % TODO: 1 .also return epochwise zscored neural firing matching
    %       2. return timeseries of smoothed firing rate
    running_times = Spk.timeBinMidPoints(Spk.sessionTypePerBin == 1);
    [behavior, thrown_out_times] = table.behavior.lookup(Option.animal, ...
                                                         running_times);
    Components = analysis.timeVarying_v2(Patterns, Option, Spk);
    Components = plots.temporal.correlateSpectral(Components, Events, Option);
    Components = plots.temporal.correlateBehavior(Components, Events, Option);
end

if Option.analysis.checks
    % Plots regarding the raw and processed data (and sometimes
    % relation to processed Patterns struct)
    % TODO: Think about splitting this into checks involving
    %        versus not involving the Patterns struct
    plots.runChecks(Events, Spk, Patterns, Option);
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
        warning("no existing DetailedRunsSummary table")
        DetailedRunsSummary = [];
    end

    % Identifying information about this options set and date of run
    hash      = DataHash(Option);
    hash      = hash(1:7); % Take the first 7 letters of the hash
    hash      = string(hash);
    timestamp = string(date());

    % Determine information to add to table
    Optim=params.getOptimizationParams(Patterns,Events,Option);
    Optimtable=struct2table(Optim);
    Optimtable.timestamp = timestamp;
    Optimtable.hash = hash;
    Optimtable.numWindowsCut = numWindowsCut;
    Optimtable.cutoffs = Events.cutoffs;

    % Create a table of options and table of patterns
    Optiontable  = struct2table(Option, 'AsArray', true);
    Patterntable = query.getPatternTable(Patterns);

    % Combine option columns with hash and date
    tablerow = [Optiontable, Optimtable];
    % Combine those with all rows of the pattern table
    tablecontent = util.table.flexibleColumnCat(Patterntable, ...
                        repmat(tablerow, height(Patterntable), 1));

    %% Check and Hash
    if ~isempty(RunsSummary)  &&any(contains(RunsSummary.hash, hash))
        RunsSummary(contains(RunsSummary.hash, hash), :) = []; % Delete any rows that contain the current hash
        DetailedRunsSummary(contains(DetailedRunsSummary.hash, hash), :) = [];
        disp("already computed before, rehashing to the same location");
        % New options:    Append row
    else
        disp("new results  --not in existing table")
    end

    % --------------------------------------
    % Append the new results and posrpocess
    % --------------------------------------
    if istable(DetailedRunsSummary)
        old_height = height(DetailedRunsSummary);
    else
        old_height = 0;
    end
    if  old_height ~= 0  && ...
        size(DetailedRunsSummary,2) ~= size(tablecontent,2)
        DetailedRunsSummary = table.addNewColumn(DetailedRunsSummary, tablecontent);
        RunsSummary         = table.addNewColumn(RunsSummary, tablerow);
    else
        DetailedRunsSummary = [DetailedRunsSummary; tablecontent];
        RunsSummary         = [RunsSummary; tablerow];
    end
    assert(height(DetailedRunsSummary) > old_height, "appending failed!");
    DetailedRunsSummary = table.postprocessSummaryTable(DetailedRunsSummary);
    RunsSummary         = table.postprocessSummaryTable(RunsSummary);

    %% ------------- Save ----------------------------
    % save the tables
    save("RunsSummary", "RunsSummary",'-v7.3');
    save("DetailedRunsSummary", "DetailedRunsSummary", '-v7.3');
    % save the results
    saveVars = {'Option'};
    if exist('Patterns','var')
        saveVars = [saveVars, {'Patterns', 'Patterns_overall'}];
    end
    if exist('Components', 'var')
        saveVars = [saveVars, {'Components', 'Components_overall'}];
    end
    thisFile = fullfile(codedefine, "hash", hash);
    disp("Saving ...");
    tic; save(thisFile, saveVars{:},'-v7.3');
    disp("... " + toc + " seconds");
    % link most recent state
    recencyName = Option.animal + "_" + replace(Option.generateH," ", "") + ...
                    "_mostRecentState";
    recencyFile = fullfile(codedefine, recencyName);
    system(['ln -sf', thisFile, ' ', recentFile]);
    % save raw?
    if Option.saveRaw
        disp("Saving raw...");
        tic; save(thisFile, "Events", "Spk",'-v7.3', '-append');
        disp("... " + toc + " seconds");
    end
    
end
