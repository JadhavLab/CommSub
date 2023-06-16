clear
%% Paths
if ~exist('pathSet','var') || pathSet == 0 % only set path if unset
    addpath(genpath(pwd)) % all folders in utils added, including semedo code
    if ispc
        paths = " C:\Users\BrainMaker\commsubspace\SingleDayExpt";
    elseif ismac
        paths(1) = "/Volumes/sharespace-commsub/data";
        paths(2) = "~/Data/commsubspace";
    else
        paths = "";
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
animal_list = ["JS21","ZT2","ER1","JS14","JS13","JS17"];
%animal_list = ["JS13","JS17"];
animals_behavior = cell(1,numel(animal_list));
animal_critical_behaviors = cell(1, numel(animal_list));
pickup_where_left_off = true;
partitioned = false


for iAnimal = progress(1:numel(animal_list),'Title','Animal')
    Option.animal = animal_list(iAnimal);
    if pickup_where_left_off == false
        continue
    end
    
    if ~partitioned % use all the HPC-PFC population
        ComponentAnalysis
        animals_behavior{iAnimal} = animal_behavior;
        animal_critical_behaviors{iAnimal} = critical_behaviors;
    else
        
    end
    %catch MatlabException
        %warning('Failed to run %s %s', Option.animal, Option.generateH);
    %end
end 
