function Default = defaults(name)
% function Default = defaults(name)
%
% Returns the default structure encoding some basic options
%
% This script causes default values to change across scripts calling it
%
% 
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

% used to determine how large the individual bins are in spikeCountMatrix
Default.spikeBinSize                  = 0.050;           % 30 milliseconds, prviously 150 ms TIME BIN FOR EACH SAMPLE OF TRRIAL

% ----- Trial and trial binning -----
Default.winSize                       = [0, 0.300];   % size of the
                                                         % full-wave trial
                                                         % window -- OVERAL
                                                         % LTRIAL WINDOW
Default.positiveDerivativeCheck       = [Default.winSize(1), Default.winSize(2)/2]; % period of window with enforced positive derivative

% used to be [-0.15, 0.15]
Default.equalWindowsAcrossPatterns    = true;            % whether all three patterns have the same #windows
Default.quantileToMakeWindows         = 0.85;

Default.thetadelta_outlierQuantile    = [0, 0.995];     % quantile to remove outliers from theta and delta

% ---- TRIAL WINDOW -----
% how long -- this will interpolate all to be timesPerTrial long
Default.spikeShiftSize                = 0.010; % used to default Default.timesPerTrial - size of actual trial sample
Default.timesPerTrial                 = ceil((sum(Default.winSize))/Default.spikeShiftSize); % 
assert(Default.timesPerTrial > 1);

% About brain areas
Default.sourceArea                    = "CA1"; % only when there are
Default.waysOfPartitions              = 2;     

% ---  Controls ---- 
Default.singleControl                 = false; % whether to use just one control column
Default.oldControlBehavior            = false; % whether to call genratePatternShuffle for control (Which creates shuffled Hc. if false, Hc=H -- which needs to be true for finding low periods of activity)
Default.lowerControl                  = true;  % whether to use the lower control (Hc) or the upper control (Hc2)

% ----- Preprocessing -----
Default.binsToMatchFR                 = 20;
Default.preProcess_FilterLowFR        = true;  % whether to filter out low FR cells
Default.preProcess_matchingDiscreteFR = true;  % whether to match discrete FR
Default.preProcess_gaussianFilter     = 5;     % gaussian that covers this many samples
Default.preProcess_zscore             = true;  % whether to zscore the data

% ----- Data resampling -----
% Data resamlping
Default.numPartition                  = 50; % subsamples of target/source split

% Time-varying
Default.dimCompAnalysis               = 5; % number of components to use in dimensionality reduction
Default.stablePerf                    = 0.9; % performance threshold for stable periods

% ----- Analysis -----
% What to run?
Default.generateH                       = join(["fromSpectra","fromRipTimes"], "  ");
Default.analysis.run_selected_genH      = false; % 📜 this
Default.analysis.checks                 = true;  % whether to run checks: FR, singleEvents, H
Default.analysis.rankRegress            = true;  % rank regression
Default.analysis.factorAnalysis         = false; % took too long -- tend to not run this
Default.analysis.singleNeuronPrediction = false; % single neuron prediction in RRR
Default.analysis.timeVarying            = false; % rewrite
Default.analysis.cca                    = false; % kim 2022
Default.analysis.JPECC                  = false; % joint peri-event cannonical correlation
Default.analysis.reverse                = false; % commsubspace high -> network patterns
% Default.generateH = "fromWpli " + " fromRipTimes";
% Default.generateH = "fromCoherence "+" fromRipTimes";
% Default.generateH = "fromSpectra "+" fromRipTimes";

% ---- Analysis parameters -----
Default.rankRegress.cvnum = 10; % number of cross-validation folds
Default.jpecc.cvnum = 4; % number of cross-validation folds

Default.save    = true; % whether to save the results
Default.saveRaw = false; % whether to save the raw data


%% --------- OPTIONS SPECIFIC TO SCRIPT TYPE ---------------------------------
%  one could place options that are specific to a script type here
if     name == "TheScript"
elseif name == "RRR_dPCA"
elseif name == "MethodTest"
end

Default = option.postamble(Default);
