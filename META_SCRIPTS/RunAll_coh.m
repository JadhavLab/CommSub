clear
commsubspaceToPath
% RunAll.m
% Runs TheScript.m for all animals and all methods
%
% Assumptions:
%   - You are in the parent ComSub directory
%   - You have added all the subdirectories to the path
%           e.g. addpath(genpath(parent_directory))

%% Script parameters
Option = option.defaults();
Option.tableAppend = "_coh";

animal_list = [...
    "JS21",...
    "JS15",...
    "JS14",...
    "JS13",...
    "JS17",...
    "ER1",...
    "ZT2"...
];

h_methods = [  ...
    ..."fromSpectra  fromRipTimes"   ...SPECTRAL POWER
    "fromCoherence  fromRipTimes"... COHERENCE
    ..."fromWpli  fromRipTimes", ...    WPLI
];

%  Load previous progress
progress_file = fullfile(hashdefine(), "last_run_coh.mat");
if exist(progress_file, "file")
    try
        load(progress_file, "last_run");
    catch
        last_run = [];
    end
else
    last_run = [];
end

%% Print what we're doing
disp(" ----------- RunAll ----------------------")
disp("Running " + numel(animal_list) + " animals");
disp("with " + numel(h_methods) + " methods");
disp("and Option struct ")
disp(rmfield(Option, {'animal', 'generateH'}));
disp("and analysis struct ")
disp(Option.analysis);
disp("Last run was " + last_run);
disp(" -----------------------------------------")
disp("Press any key to continue");
pause
first = false;

dopar = false;
if dopar
    jobs = [];
end

%% Run
pickup_where_left_off = false; % false to look for last processed, true to default process all (loop sets a false to true when finds correct last animal/genH)
[cntAn, cntH] = deal(0);
for genH_= progress(h_methods,'Title','genH method'); cntH = cntH + 1;
for iAnimal = progress(1:numel(animal_list),'Title','Animal'); cntAn = cntAn + 1;
        Option.animal = animal_list(iAnimal);
        Option.generateH = genH_;
        if ~isequal([Option.animal, Option.generateH], last_run)
           pickup_where_left_off = true;
        else
           disp("Skipping previous run");
        end
        if pickup_where_left_off == false
            continue
        end
        disp(newline + "-------------------------------");
        disp("Running " + Option.animal + " " + Option.generateH);
        disp("-------------------------------");
        if ~first
            pause
            disp("Press any key to continue");
            first = true;
        end
        if dopar
            jobs(cntH, cntAn) = batch(TheScript, ...
            'Workspace', {Option}, 'CurrentFolder', pwd);
        else
            %try
            diary(figuredefine("logs", replace(strjoin([Option.animal, Option.generateH], "_"), " ", "") + ".log"));
            TheScript
            diary off
            %catch MatlabException
                %warning('Failed to run %s %s', Option.animal, Option.generateH);
            %end
            last_run = [Option.animal, Option.generateH];
            save(progress_file, "last_run");
        end
end % genH
end % animal
