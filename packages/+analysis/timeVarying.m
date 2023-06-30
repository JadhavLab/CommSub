function Components = timeVarying(Patterns, Option, r)
% Components = timeVarying(Patterns, Option)
%
% This function takes the patterns and performs the time varying analysis
% on them. It returns a struct with the critical components for each
% behavior and the subspaces for each behavior.
%
% Inputs:
%   Patterns: struct with the patterns for each partition
%   Option: struct with the options for the analysis
%
% Outputs:
%   Components: struct with the critical components for each behavior and
%   the subspaces for each behavior


disp("Running time varying analysis")
tic

const = Option.shortcut;

if ~isfield(Patterns, 'rankRegress') || isempty(Patterns(1).rankRegress)
    error("Patterns must have the rankRegress field")
end

%% Take the behavior table
running_spikeTimes = r.timeBinMidPoints(r.sessionTypePerBin == 1);
[animal_behavior, throwout_times] = table.behavior.lookup(Option.animal, ...
    running_spikeTimes);
[animal_behavior,unique_times] = behaviors.addBehToTable(animal_behavior);

clear Components
Components= struct(...
    "compAnalysisByBeh",[], ...
    "allSubSpaceComponents",[]...
);
Components = repmat(Components, ...
                    [Option.numPartition, Option.waysOfPartitions]);

for p = 1:Option.numPartition
    for j = [HPC, PFC]

        target = const.areanames(j);
        
        B_ = cell(1,Option.nPatternAndControl);
        for i = 1:Option.nPatternAndControl
            B_{i} = Patterns(p,j,i).rankRegress.B_;
        end
        
        [total_subspaces, cell_subspaces] = components.subspaceSimilarity...
            (Option.dimCompAnalysis, B_, r.spikeRateMatrix, r.celllookup, Option.sourceArea, target,...
            'source_index', Patterns(p,j,i).index_source,...
            'target_index', Patterns(p,j,i).index_target);
        
        [behavior_running_subspaces] = components.animalBehComponents...
            (unique_times, throwout_times, total_subspaces,sessionTypePerBin);
 
        [critical_behaviors, unifiedTime, critical_components, critical_times] = ...
            components.makeBehaviorStructs(Option.dimCompAnalysis, cellOfWindows, ...
            animal_behavior, behavior_running_subspaces, false);

        Components(p,j).compAnalysisByBeh = critical_behaviors;
        Components(p,j).allSubSpaceComponents = behavior_running_subspaces;
        for k = 1:4
            Components(p,j).compAnalysisByBeh(k).critical_components = critical_components{k};
            Components(p,j).compAnalysisByBeh(k).critical_times = critical_times{k};
            Components(p,j).compAnalysisByBeh(k).unifiedTime = unifiedTime;
        end
    end
end

disp("Finished time varying analysis in " + string(toc) + " seconds")


