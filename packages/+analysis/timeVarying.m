function Components = timeVarying(Patterns, Option)

%% Take the behavior table
running_spikeTimes = timeBinMidPoints(sessionTypePerBin == 1);
[animal_behavior, throwout_times] = table.behavior.lookup(animal, running_spikeTimes);
[animal_behavior,unique_times] = behaviors.addBehToTable(animal_behavior);

clear Components
Components= struct("compAnalysisByBeh",[], "allSubSpaceComponents",[]);
Components = repmat(Components, ...
    [Option.numPartition, Option.waysOfPartitions]);

for p = 1:Option.numPartition
    for j = [HPC, PFC]
        if j == HPC
            target = "CA1";
        else
            target = "PFC";
        end
        
        B_ = cell(1,6);
        for i = 1:nPatternAndControl
            B_{i} = Patterns(p,j,i).rankRegress.B_;
        end
        
        [total_subspaces, cell_subspaces] = components.subspaceSimilarity...
            (Option.dimCompAnalysis, B_, spikeRateMatrix, celllookup, Option.sourceArea, target,...
            'source_index', Patterns(p,j,i).index_source,...
            'target_index', Patterns(p,j,i).index_target);
        
        [ behavior_running_subspaces] = components.animalBehComponents...
            (unique_times, throwout_times, total_subspaces,sessionTypePerBin);
        % also add the subspace strength to the table
 
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


