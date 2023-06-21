function Default = defaults(name)
% function Default = defaults(name)
%
% Returns the default structure encoding some basic options
%
% Now changes in default values automatically occur across scripts
%
% Important variables we iterate over:
%   - Default.generateH = "fromCoherence "+" fromRipTimes" |
%                         "fromSpectra "+" fromRipTimes"
%   - Default.animal = "ZT2" | "JS13" | "JS14" | "JS21" | "ER1" | "JS17"
%   componet analysis: animal_list = ["JS21","ZT2","ER1","JS14","JS13","JS17"];
%
%  Less important iteration variables:
%  - Default.sourceArea = "CA1" | "PFC" % not enough PFC cells usually for
%                                       % source to be PFC


if nargin == 0
    name = "TheScript";
end


%% Script parameters
Default = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----- Key properties -----
Default.animal                        = "ZT2";
Default.samplingRate                  = [] ;             % For spikes.getSpikeTrain, nan if not given
% ----- Trial and trial binning -----
Default.winSize                       = [-0.15, 0.15];   % size of the
                                                         % full-wave trial
                                                         % window -- OVERAL
                                                         % LTRIAL WINDOW
Default.spikeBinSize                  = 0.050;           % 30 milliseconds, prviously 150 ms TIME BIN FOR EACH SAMPLE OF TRRIAL
Default.spikeShiftSize                = 0.010;           % 30 milliseconds, prviously 150 ms HOW MUCH TO SHIFT TO CRERAT EACH SAMLPE
Default.timesPerTrial                 = ceil(range(Default.winSize)/Default.spikeShiftSize);
assert(Default.timesPerTrial > 1);
Default.equalWindowsAcrossPatterns    = true;            % whether all three patterns have the same #windows
Default.quantileToMakeWindows         = 0.85;

% About brain areas
Default.sourceArea                    = "CA1";           % only when there are
Default.waysOfPartitions              = 2;

% Controls
Default.singleControl                 = false;           % whether to use just one control column
Default.oldControlBehavior            = false;
Default.lowerControl                  = true;
Default.binsToMatchFR                 = 20;
Default.preProcess_FilterLowFR        = true;
Default.preProcess_matchingDiscreteFR = true;

% Data resamlping
Default.numPartition                  = 50;              % ways to split source and target

% Time-varying
Default.dimCompAnalysis               = 5;
Default.stablePerf                    = 0.9;

% What to run?
Default.generateH                = join(["fromSpectra","fromRipTimes"], "  ");
Default.analysis.run_selected_genH      = false;
Default.analysis.frChecks               = true;
Default.analysis.rankRegress            = true;
Default.analysis.factorAnalysis         = false;
Default.analysis.timeVarying            = false;
Default.analysis.singleNeuronPrediction = false;
% Default.generateH = "fromWpli " + " fromRipTimes";
% Default.generateH = "fromCoherence "+" fromRipTimes";
% Default.generateH = "fromSpectra "+" fromRipTimes";
%
Default.rankRegress.cvnum = 10; % number of cross-validation folds
Default.jpecc.cvnum = 4; % number of cross-validation folds

Default.save    = true; % whether to save the results
Default.saveRaw = false; % whether to save the raw data


%% --------- OPTIONS SPECIFIC TO SCRIPT TYPE ---------------------------------
if     name == "TheScript"
elseif name == "RRR_dPCA"
elseif name == "MethodTest"
end

Default = option.postamble(Default);
