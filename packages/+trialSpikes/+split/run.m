function out = run(r, Option)
% trialSpikes.split.run
%
% Input:
%   r: trialSpikes object
%   Option: struct
%       .binsToMatchFR: number of bins to match FR
%       .sourceArea: 'hpc' or 'pfc'
%
% Output:
%   out: struct
%       .source: source spikes
%       .target: target spikes
%       .nSource: number of source cells
%       .nTarget: number of target cells
%       .source_index: index of source cells
%       .target_index: index of target cells

    out.sourceArea = Option.sourceArea;

    if strcmp(Option.sourceArea, "CA1")
        source = r.hpc.X;
        other_target = r.pfc.X;
        out.directionality = ["hpc-hpc", "hpc-pfc"];
    elseif strcmp(Option.sourceArea, "PFC")
        source = r.pfc.X;
        other_target = r.hpc.X;
        out.directionality = ["pfc-pfc", "pfc-hpc"];
    else
        error("Invalid source area");
    end

    [source, target, nSource, nTarget, source_index, target_index] ...
        = trialSpikes.split.twoWayFRMatch(r.hpc.X, r.pfc.X, ...
                      r.hpc.FR, r.pfc.FR, Option.binsToMatchFR);
    out.source = source;
    out.target = target;
    out.nSource = nSource;
    out.nTarget = nTarget;
    out.source_index = source_index;
    out.target_index = target_index;
    out.source_area = Option.sourceArea;

end

