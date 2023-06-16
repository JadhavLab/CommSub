clear
%% Paths
if ~exist('pathSet','var') || pathSet == 0 % only set path if unset
    addpath(genpath(pwd)) % all folders in utils added, including semedo code
    if ispc
        paths = " C:\Users\BrainMaker\commsubspace\SingleDayExpt";
    elseif ismac
        paths(1) = "/Volumes/sharespace-commsub/data";
        paths(2) = "~/Data/commsubspace";
    end

    arrayfun(@(path) addpath(genpath(path)), paths);
    pathSet = 1;
end
% %%
% FIELD = [];
% numToRun = numel(FIELD);
%% Script parameters
% -----------------------------------------------------------------
Option = struct();

% Option.animal =

% Option.generateFromRipTimes = true; 
% Option.generateH = "fromFilteredEEG "+" fromRipTimes";
% Option.samplingRate  = "empty" ;         % For spikes.getSpikeTrain
% Option.spikeBinSize  = 0.1;          % 100 milliseconds
% Option.timesPerTrial = 10;         % 10 times per trial
% Option.winSize       = {[-0.15, 0.15]};              % size of the window
% Option.sourceArea    = "PFC";
% Option.equalWindowsAcrossPatterns = true;    % whether all three patterns have the same #windows
% Option.singleControl = true;                 % whether to use just one control column
% Option.numPartition = 10;                    % ways to split source and target
% usingSingleprediction = true;
% 
% winSize = Option.winSize{1};

%animal_list = ["JS21","JS15","JS14","JS13", "JS17"];
animal_list = ["ER1"];

pickup_where_left_off = true;
for genH_= progress(["fromWpli  fromRipTime",],'Title','genH method')
for iAnimal = progress(1:numel(animal_list),'Title','Animal')
    Option.animal = animal_list(iAnimal);
    Option.generateH = genH_;
    %if isequal([Option.animal, Option.generateH], ...
    %        ["JS13", "fromCoherence  fromRipTime"])
    %    pickup_where_left_off = true;
    %else
    %    disp("Skipping previous run");
    %end
    if pickup_where_left_off == false
        continue
    end
    %try
    TheScript
    %catch MatlabException
        %warning('Failed to run %s %s', Option.animal, Option.generateH);
    %end
end 
end
