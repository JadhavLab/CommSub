function Out = timeVarying(Patterns, Option, r)
% Components = timeVarying(Patterns, Option)
%
% This function takes the patterns and performs the time varying analysis
% on them. It returns a struct with the critical components for each
% behavior and the subspaces for each behavior.
%
% TODO: 
% 1. How much each cell matches the subspace
% 2. How much all match
%    - old RRR method
%    - new RRR method
%    - CCA method
% 3. Relationship to various behavioral variables
%   - speed
%   - acceleration
%   - choice point (out in)
%   - reward location (out in)
%
% Inputs:
% -----
%   Patterns: struct with the patterns for each partition
%   Option: struct with the options for the analysis
%
% Outputs:
%   Components: struct with the critical components for each behavior and
%   the subspaces for each behavior

% ip = inputParser();
% ip.addParameter('animal_behavior', []);
% % ip.addParameter('unique_times', []);
% % ip.addParameter('throwout_times', []);
% ip.parse(varargin{:});
% Option = ip.Results;

disp("Running time varying analysis")
tic
Const = option.constants();

if ~isfield(Patterns, 'rankRegress') || isempty(Patterns(1).rankRegress)
    error("Patterns must have the rankRegress field")
end

% %% Take the behavior table
% if isempty(Option.animal_behavior)
%     running_spikeTimes = r.timeBinMidPoints(r.sessionTypePerBin == 1);
%     [animal_behavior, throwout_times] = table.behavior.lookup(Option.animal, ...
%         running_spikeTimes);
%     [animal_behavior, unique_times] = behaviors.addBehToTable(animal_behavior);
% else
%     animal_behavior = Option.animal_behavior;
% end

clear Components
Out= struct(...
    'rankRegress', [],...
    'target', []...
    );
Out = repmat(Out, [Option.numPartition, Option.waysOfPartitions]);
% Components = repmat(Components, ...
%                     [Option.numPartition, Option.waysOfPartitions]);

for p = 1:min(Option.numPartition, size(Patterns, 1))
    for j = [Const.HPC, Const.PFC]

        O = Out(p, j);
        P = Patterns(p, j);

        O.target = Const.areanames(j);

        % old method
        % ----------------
        % old = analysis.temporal.oldMatching(p, Option, r, target, animal_behavior,...
        %                               unique_times, throwout_times);

        % new method
        % ----------------
        for method = {'prod', 'concat'}
            O.rankRegress.(method{1}) = analysis.temporal.match(P, Option, r,...
                        'method', method{1},...
                        'component_method', 'rankRegress'...
                        );
        end

        Out(p, j) = O;
    end
end

disp("Finished time varying analysis in " + string(toc) + " seconds")
