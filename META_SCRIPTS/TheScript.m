
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
% Documentation
% Events is a struct with fields:
% - .times : array of times of events
% - .H     : Event Matrix, T x 3, and each column are theta, delta, ripple
% - .Hvals : Event Matrix, T x 3, values without nans
% - .Hnanlocs : Event Matrix, T x 3, logicals of nans

%%%%%%%%%%%%%%%% CUT WINDOWS WITH EVENT MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[cellOfWindows, cutoffs] = windows.ThetaDeltaRipple(Events, Option);
% Documentation
% -  cellOfWindows: 1 x nPatterns cell array of windows
% -  cutoffs: nPatterns x 1 vector of cutoffs

%%%%%%%%%%%%%%%% ACQUIRE SPIKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Getting spikes
[timeBinStartEnd, timeBinMidPoints, ~, spikeCountMatrix, spikeRateMatrix, ...
    areaPerNeuron, cell_index, sessionTypePerBin] = spikes.getSpikeTrain(Option.animal, ...
    Option.spikeBinSize,  Option.samplingRate);

% filter the neurons whose firing rate is lower than specified threshold
if Option.preProcess_FilterLowFR 
    disp("Filtering low FR neurons")
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
r.celllookup = cellInfo.getCellIdentities(Option.animal, cell_index,...
                                          areaPerNeuron);
r.avgFR         = avgFR;
r.areaPerNeuron = areaPerNeuron;
r.pfc.nNeurons = nPFCneurons;
r.hpc.nNeurons = nHPCneurons;
r.windowInfo.cellOfWindows = cellOfWindows;
r.windowInfo.nWindows      = cellfun(@(x) size(x, 1), cellOfWindows);
r.nPattern = 3;
r.nControl = Option.nPatternAndControl - r.nPattern;
r.timeBinMidPoints = timeBinMidPoints;
r.sessionTypePerBin = sessionTypePerBin;
frChecks(r, "appendFigTitle", r.animal);


%%%%%%%%%%%%%%%% SETUP PARTITIONS AND RESULT STRUCTURES %%%%%%%%%%%%%%%%%%
Patterns = trialSpikes.partitionAndInitialize(r, Option);

%%%%%%%%%%%%%%%% ANALYSIS SECTION    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Option.analysis.rankRegress
    % TODO: 
    % 1. fix Option.rankregress => Option.rankRegress
    % 2. most rankRegress.B_ are empty
    Patterns = analysis.rankRegress(Patterns, Option);
end

if Option.analysis.timeVarying
    % ISSUE: hits a bug on line 4
    % TODO: 1 .also return epochwise zscored neural firing matching
    %       2. return timeseries of smoothed firing rate
    Components = analysis.timeVarying(Patterns, Option);
end

if Option.analysis.factorAnalysis
    Patterns = analysis.factorAnalysis(Patterns, Option);
end

%%%%%%%%%%%% Assign all of the RAW relevent structures %%%%%%%%%%%%
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

%%%%%%%%%%%%%%%% CREATE TABLE AND SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare table!
path2 = "C:\Users\BrainMaker\MATLAB Drive\Shared";
cd(path2)

%% Megatable: All of the information for each partition run above
if exist("Megatable.mat", 'file')
    load("Megatable.mat");
else
    TABLE = [];
end

%% T: Summary information (not partition-wise) for this entire set of runs
if exist("T.mat", 'file')
    load("T.mat");
else
    warning("no exisiting real table")
    T = [];
end

params.optimizeOptionsScript;
%Option = rmfield(Option, "animal");
Optiontable  = struct2table(Option, 'AsArray', true);
Patterntable = query.getPatternTable(Patterns);

factorAnalysisChecksum = sum(Patterntable.singularWarning) ...
                         > 0.7*size(Patterntable,1);
if factorAnalysisChecksum
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
if size(T,2) ~= size(tablecontent,2)
    T = table.addNewColumn(T, tablecontent);
    TABLE    = table.addNewColumn(TABLE, tablerow);
else
    T = [T; tablecontent];
    TABLE = [TABLE; tablerow];
end
assert(height(T) > old_height, "not appending");

%% ------------- Save ----------------------------
save("Megatable", "TABLE",'-v7.3');
save("T", "T", '-v7.3');
% writetable(TABLE,'~/Megatable.xls');
%save(fullfile(datadefine, "hash_forTesting", hash), ...
 %   "Patterns", "Option", "animal_behavior",'Behaviors',"celllookup",'-v7.3')

save(fullfile(datadefine, "hash_forTesting", hash), ...
    "Option", "animal_behavior",'Behaviors',"celllookup",'-v7.3')
% most recent state
save(fullfile(datadefine, 'mostRecentState'))

