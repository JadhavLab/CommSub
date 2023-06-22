function [s_hpc, s_pfc, t_hpc, t_pfc, s_hpc_index, s_pfc_index, t_hpc_index,...
          t_pfc_index, directionality] = legacyRun(r, Option)
% Partition (subsample) the available cells into source and target ensuring
% that the firing rate distribution of the source and target are matched.
%
% Inputs
% ------
% r : struct
%   r.hpc : struct
%       r.hpc.FR : double
%           overall firing rate of hpc neurons
%       r.hpc.X : double
%           firing rate matrix of hpc neurons
%   r.pfc : struct
%       r.pfc.FR : double
%           overall firing rate of pfc neurons
%       r.pfc.X : double
%           firing rate matrix of pfc neurons
% Option : struct
%   Option.waysOfPartitions : int
%       number of ways to split the data
%   Option.nPatternAndControl : int
%       number of neurons to use for pattern and control
%   Option.sourceArea : string
%       source area for the split
%   Option.targetArea : string
%       target area for the split
%
% Outputs
% -------
% s_hpc : double
% s_pfc : double
% t_hpc : double
% t_pfc : double
% s_hpc_index : double
% s_pfc_index : double
% t_hpc_index : double
% t_pfc_index : double

s_hpc = nan; s_pfc = nan; s_hpc_index = nan; s_pfc_index = nan;

switch Option.waysOfPartitions

    case 4 % 4 WAY SPLIT
        [s_hpc, s_pfc, t_hpc, t_pfc, s_hpc_index, s_pfc_index, t_hpc_index,t_pfc_index ] = ...
            trialSpikes.split.fourWay(r.pfc.FR, r.hpc.FR, r.hpc.X, r.pfc.X, Option.binsToMatchFR);
        directionality = ["hpc-hpc","hpc-pfc", "hpc-hpc", "pfc-hpc"];
        return;
    case 3 % 3 WAY SPLIT (three-way, where source split in two)
        if Option.targetArea == "CA1"
            [s_hpc, s_pfc, t_ca1, s_hpc_index, s_pfc_index, t_ca1_index]  = ...
            trialSpikes.split.threeWay(r.hpc.FR, r.pfc.FR, r.hpc.X, r.pfc.X, ...
            Option.binsToMatchFR);
        else
            [s_pfc, s_hpc, t_pfc, s_pfc_index, s_hpc_index, t_pfc_index] = ...
            trialSpikes.split.threeWay(r.pfc.FR, r.hpc.FR, r.pfc.X, r.hpc.X,...
            Option.binsToMatchFR);
        end
        if targetArea == "CA1"
            directionality = ["hpc-hpc","pfc-hpc"];
        else
            directionality = ["pfc-hpc", "pfc-pfc"];
        end

    case 2 % TWO WAY SPLIT (three-way, where target split in two)
        if Option.sourceArea == "CA1"
            directionality = ["hpc-hpc","hpc-pfc"];
        else
            directionality = ["pfc-hpc","pfc-pfc"];
        end
        pure_random_select = ~Option.preProcess_matchingDiscreteFR;
        if pure_random_select
            % (DEPRICATED) : not really in use as much
            % randomly select neurons based on number of neurons in the
            % respective regions and split source and target
            [X_source,X_target, nSource, nTarget,...
                index_source, index_target] = ...
                trialSpikes.split.sourceTarget(2, ...
                Option.nPatternAndControl, r.hpc.X, r.pfc.X, ...
                'specifiedSource', Option.sourceArea);
        else
            % (PRIMARY METHOD)
            % match the firing rate distribution
            if lower(Option.sourceArea) == "ca1"
                [s_hpc, target, nSource, nTarget, s_hpc_index, index_target] ...
                    = trialSpikes.split.twoWayFRMatch(r.hpc.X, ...
                        r.pfc.X, r.hpc.FR, r.pfc.FR, Option.binsToMatchFR);
                t_hpc = target(1, :, :);
                t_pfc = target(2, :, :);
                t_hpc_index = index_target(1, :, :);
                t_pfc_index = index_target(2, :, :);
            elseif lower(Option.sourceArea) == "pfc"
                [s_pfc,target, nSource, nTarget,...
                    s_pfc, index_target] = ...
                    trialSpikes.split.twoWayFRMatch(r.hpc.X, r.pfc.X, ...
                                r.hpc.FR, r.pfc.FR, Option.binsToMatchFR);

                t_pfc = target(1,:,:);
                t_hpc = target(2,:,:);
                t_pfc_index = index_target(1,:,:);
                t_hpc_index = index_target(2,:,:);
            else
                error('Option not implemented')
            end
        end
    otherwise
        error('Option not implemented')
end
