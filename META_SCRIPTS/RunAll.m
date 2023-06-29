clear
% RunAll.m
% Runs TheScript.m for all animals and all methods
%
% Assumptions:
%   - You are in the parent ComSub directory
%   - You have added all the subdirectories to the path
%           e.g. addpath(genpath(parent_directory))

%% Script parameters
Option = option.defaults();
% -----------------------------------------------------------------
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
% winSize = Option.winSize{1};

animal_list = [...
                "JS21",...
                "JS15",...
                "JS14",...
                "JS13",...
                "JS17",...
                "ER1"...
                ..."ZT2"...
];
h_methods = [  ...
             "fromSpectra  fromRipTimes"   ...SPECTRAL POWER
             "fromCoherence  fromRipTimes"... COHERENCE
             "fromWpli  fromRipTimes", ...    WPLI
            ];

%  Load previous progress
progress_file = fullfile(codedefine(), "hash", "last_run");
if exist(progress_file, "file") == 2
    last_run = load(progress_file, "last_run");
else
    last_run = [];
end

%% Print what we're doing
disp(" ----------- RunAll ----------------------")
disp("Running " + numel(animal_list) + " animals");
disp("with " + numel(h_methods) + " methods");
disp("and Option struct ")
disp(Option);
disp("and analysis struct ")
disp(Option.analysis);
disp("Last run was " + last_run);
disp(" -----------------------------------------")
disp("Press any key to continue");
pause

%% Run
pickup_where_left_off = true;
for genH_= progress(h_methods,'Title','genH method')
    for iAnimal = progress(1:numel(animal_list),'Title','Animal')
        Option.animal = animal_list(iAnimal);
        Option.generateH = genH_;
        disp("Running " + Option.animal + " " + Option.generateH);
        if isequal([Option.animal, Option.generateH], last_run)
           pickup_where_left_off = true;
        else
           disp("Skipping previous run");
        end
        if pickup_where_left_off == false
            continue
        end
        %try
        TheScript
        %catch MatlabException
            %warning('Failed to run %s %s', Option.animal, Option.generateH);
        %end
        last_run = [Option.animal, Option.generateH];
        save(progress_file, "last_run");
    end 
end
