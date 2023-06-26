function [Components] = oldMatching(Option, Patterns, r, target, ...
    animal_behavior, unique_times, throwout_times)
    % This function is used to match the components of the behavior
    % analysis to the components of the pattern analysis. Used for the subpsace
    % temporal analysis in Ziyi's thesis.
    % Input:
    %   Option: the option struct
    %   Patterns: the pattern struct
    %   r: the struct of the data
    %   const: the struct of the constants
    
    B_ = cell(1,Option.nPatternAndControl);
    for i = 1:Option.nPatternAndControl
        B_{i} = Patterns(i).rankRegress.B_;
    end
    
    [total_subspaces, ~] = ...
    components.subspaceSimilarity(Option.dimCompAnalysis, B_, ...
    r.spikeRateMatrix, r.celllookup, Option.sourceArea, target,...
        'source_index', Patterns(i).index_source,...
        'target_index', Patterns(i).index_target);
    
    [behavior_running_subspaces] =...
    components.animalBehComponents(unique_times, throwout_times,...
    total_subspaces,r.sessionTypePerBin);

    [critical_behaviors, unifiedTime, critical_components, critical_times] ...
    = components.makeBehaviorStructs(Option.dimCompAnalysis, r.cellOfWindows, ...
    animal_behavior, behavior_running_subspaces, false);

    Components.compAnalysisByBeh = critical_behaviors;
    Components.allSubSpaceComponents = behavior_running_subspaces;
    for k = 1:4
        Components.compAnalysisByBeh(k).critical_components = ...
                                                 critical_components{k};
        Components.compAnalysisByBeh(k).critical_times = critical_times{k};
        Components.compAnalysisByBeh(k).unifiedTime = unifiedTime;
    end

end
