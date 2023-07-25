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
    disp("Option struct does not exist, creating one")
    addpath(genpath(codedefine()));
    Option = option.defaults(); 
else
    disp("Option struct already exists, using that")
    Option = option.setdefaults(Option);
    disp("Option struct is: ")
    disp(Option)
end
%%%%%%%%%%%%%%%% DISPLAY OUR OPTIONS TO USER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(Option.loadifexists, false) && ...
    exist(store.gethash(Option) + ".mat", 'file')
    disp("Loading from file: " + store.gethash(Option) + ".mat")
    m = matfile(store.gethash(Option) + ".mat");
    % m = matfile("bef0923.mat", "Writable", true);
    disp("Loaded variables: ")
    Events             = util.matfile.getdefault(m, 'Events', []);
    Spk                = util.matfile.getdefault(m, 'Spk', []);
    Patterns           = util.matfile.getdefault(m, 'Patterns', []);
    Patterns_overall   = util.matfile.getdefault(m, 'Patterns_overall', []);
    Components         = util.matfile.getdefault(m, 'Components', []);
    Components_overall = util.matfile.getdefault(m, 'Components_overall', []);
    Option             = util.matfile.getdefault(m, 'Option', []);
    disp("...done")
else
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
    % Getting spikes
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
        if ~isfield(Spk, 'muFR')
            Spk.muFR  = mean(Spk.spikeRateMatrix, 2);
            Spk.stdFR = std(Spk.spikeRateMatrix, 0, 2);
        end
        Spk.spikeRateMatrix  = zscore(Spk.spikeRateMatrix,  0, 2);
        Spk.spikeCountMatrix = zscore(Spk.spikeCountMatrix, 0, 2);
        Spk.avgFR = mean(Spk.spikeRateMatrix, 2);
    end
    prewindow_copy = Spk;

    % %%%%%%%%%%%%%% ACQUIRE TRIALS FROM WINDOWS + SPIKES %%%%%%%%%%%%%%%%%%%
    % RYAN bug here .. timeBinStartEnd instead of timeBinMidPoints
    disp("------------------------")
    disp("   Windowing spikes     ")
    disp("------------------------")
    Spk = trialSpikes.generate(Spk, Events, Option);

    %%%%%%%%%%%%%%%%% SETUP RAW DATA STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Structure for separated data
    %%%%%%%%%%%%%%%% SEPRATE BRAIN AREA DATASULT STRUCTURES %%%%%%%%%%%%%%%%%%
    % Separate spikesSampleMatrix/Tensor by area that neurons are in PFC and
    % neurons that in HPC
    %% Separate firing pattern into source and target
    [Spk.nSource,~,~] = size(Spk.hpc.X{1});
    [Spk.nTarget,~,~] = size(Spk.pfc.X{1});
    Spk.celllookup = cellInfo.getCellIdentities(Option.animal, Spk.cell_index,...
                                                Spk.areaPerNeuron);
    system("pushover-cli 'Finished munging data for analysis'");

    %%%%%%%%%%%%%%%% SETUP PARTITIONS AND RESULT STRUCTURES %%%%%%%%%%%%%%%%%%
    disp("------------------------")
    disp(" Subsampling partitions ")
    disp("------------------------")
    [Patterns, Patterns_overall] = trialSpikes.partitionAndInitialize(Spk, Option);
    Components = nd.initFrom(Patterns, ...
    {'index_source', 'index_target', 'directionality', 'name'});
    Components_overall = nd.initFrom(Patterns_overall, ...
    {'index_source', 'index_target', 'directionality', 'name'});
end

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
    Patterns         = analysis.rankRegress(Patterns, Option, 'verbose', true);        
    Patterns_overall = analysis.rankRegress(Patterns_overall, Option);
end

if Option.analysis.factorAnalysis
    % Factor analysis of network pattern windows of spiking activity
    % (Used to measure instrinsic dimensionality of network activity)
    Patterns = analysis.factorAnalysis(Patterns, Option);
end

if Option.analysis.cca
    running_times = Spk.timeBinMidPoints(Spk.sessionTypePerBin == 1);
    [behavior, thrown_out_times] = table.behavior.lookup(Option.animal, ...
                                                         running_times);
    % ---------------------------------------------------------------------
    % (Subspaces acquired here)
    % Patterns         = analysis.cca(Patterns, Option);
    Patterns_overall = analysis.cca(Patterns_overall, Option);
    % Copy over the components
    Components_overall =  ...
         nd.fieldSet(Components_overall, 'cca', Patterns_overall);
    % ---------------------------------------------------------------------
    % Event analysis ------------------------------------------------------
    % (append cca commsub levels during events)
    event_anal   = ... 
         analysis.cca.event_analysis(Patterns_overall, Spk, Events, Option);
    Components_overall = ... 
         nd.fieldSet(Components_overall, 'event_anal', event_anal);
    % ---------------------------------------------------------------------
    % Create table of results ---------------------------------------------
    % (create a table regarding cca versus efizz and behavior)
    efizz = load(Option.animal + "spectralBehavior.mat", "efizz");
    efizz = efizz.efizz;
    table.analyses.ccatime(Patterns_overall, efizz, Option, behavior,...
                          'behaviorColumns', {'vel', 'accel', 'lindist', 'rewarded', 'trajbound','inBoundChoiceTimes','outBoundChoiceTimes','rewardTimes'});
    % ---------------------------------------------------------------------
    % Triggered spectrogram -----------------------------------------------
    % (create compute triggered spectrograms for commsubs)
    disp("Running triggered spectrogram - run")
    close all
    triggered_spectrogram_run = ...
         analysis.cca.triggered_spectrogram(Patterns_overall, Spk, efizz,...
            'ploton', true, ... 
            'figAppend', strjoin([Option.animal,Option.genH_name], "_"), ...
            'runtype', 1);
    % disp("Running triggered spectrogram - sleep")
    % triggered_spectrogram_sleep = ...
    %      analysis.cca.triggered_spectrogram(Patterns_overall, Spk, efizz,...
    %         'ploton', true, ... 
    %         'figAppend', strjoin([Option.animal,Option.genH_name], "_"), ...
    %         'runtype', 0);
    % for i = 1:numel(Patterns_overall)
    %     Patterns_overall(i).triggered_spectrogram_run   = triggered_spectrogram_run(i,:);
    %     Patterns_overall(i).triggered_spectrogram_sleep = triggered_spectrogram_sleep(i,:);
    % end
    % ---------------------------------------------------------------------
end

if Option.analysis.timeVarying
    running_times = Spk.timeBinMidPoints(Spk.sessionTypePerBin == 1);
    [behavior, thrown_out_times] = table.behavior.lookup(Option.animal, ...
                                                         running_times);
    % How much spiking moment to moment is explained by subspace
    % Requirements: Option.analysis.cca
    % TODO: 1 .also return epochwise zscored neural firing matching
    %       2. return timeseries of smoothed firing rate
    % Component matching over time
    rrr = analysis.match_rrr(Patterns_overall, Option, Spk);
    cca = analysis.match_cca(Patterns_overall, Option, Spk);
    % Spectral matches
    Components_overall = ... 
    plots.temporal.correlateSpectral(Components_overall, Events, Option);
    Components_overall = ... thescript
    plots.temporal.correlateSpectral(Components_overall, Events, Option, 'componentMethod', 'cca');
    % Behavior matches
    Components         = plots.temporal.correlateBehavior(Components, Events, Option);
end

if Option.analysis.checks
    % Plots regarding the raw and processed data (and sometimes
    % relation to processed Patterns struct)
    % TODO: Think about splitting this into checks involving
    %        versus not involving the Patterns struct
    if strcmp(Option.animal, "JS21"); wait_state = true;
    else; wait_state = false; end
    plots.runChecks(Events, Spk, Patterns, Option, ...
                    'parallel', true, 'wait', wait_state);
end

% TODO: (1) plug in JPECC version of rankRegress here
% TODO: (2) function that outputs average response of Pattern struct per neuron

%%%%%%%%%%%%%%%% CREATE TABLE AND SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Option.save

    disp("Saving results")
    store.savetables(Events, Patterns, Option);

    saveVars = [];
    if exist('Patterns','var')
        saveVars.Patterns = Patterns;
        saveVars.Patterns_overall = Patterns_overall;
    end
    if exist('Components', 'var')
        saveVars.Components         = Components;
        saveVars.Components_overall = Components_overall;
    end
    store.savevars(Option, Events, Spk, saveVars);
end
