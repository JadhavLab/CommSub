function [Patterns, Patterns_overall] = partition(r, Option) 
% partition - Partition data into source and target areas
% 
% Parameters:
%   r - A struct containing raw and processed data
%   Option - A struct containing options for the partitioning
% 
% Returns:
%   Patterns - A struct containing the partitioned data (partitions of
%               neurons; subsampled sets of them)
%   Patterns_overall - A struct WITHOUT parititioning neurons
%
% ----------------------------
% Note:
% ----------------------------
% when the partition is three-ways, direction==1 means same target/source pair
% and direction==2 means diff target/source pair
% ----------------------------
disp('Partitioning data...')
tic
 
% Initialize the scaffold of the pattern struct
Patterns         = initPatternStruct();
Patterns_overall = initPatternStruct();

% ------------------------------
% Place paritioned data properly
% ------------------------------
for iPartition = progress(1:Option.numPartition, 'Title', 'Partitioning data')
% for iPartition = progress(1:2, 'Title', 'Partitioning data')

    % Split cells into source and target areas
    % [s_hpc, s_pfc, t_hpc, t_pfc, s_hpc_index,...
    %  s_pfc_index, t_hpc_index, t_pfc_index, directionality] = ...
    %     trialSpikes.split.legacyRun(r, Option);
    split_info = trialSpikes.split.run(r, Option);

    for i = 1:numel(r.windowInfo.cellOfWindows)
        for j = 1:numel(split_info.directionality)

            
            % Parse directionality
            directionality = split_info.directionality(j);
            sourcetarg = directionality.split('-');
            source = sourcetarg(1);
            target = sourcetarg(2);
            s_dat = split_info.source{i};
            s_ind = split_info.source_index;
            t_dat = split_info.target{j, 1, i};
            t_ind = split_info.target_index(j, :);
            
            % Assign x_source and x_target
            Patterns(iPartition,j,i).X_source  = s_dat;
            Patterns(iPartition,j,i).X_target  = t_dat;
            Patterns(iPartition,j,i).X_time    = reshape(r.trialTimes{i}',1,[]);
            % Assign index_source and index_target
            Patterns(iPartition,j,i).index_source = s_ind;
            Patterns(iPartition,j,i).index_target = t_ind;
            % Assign directionality
            Patterns(iPartition,j,i).directionality = directionality;
            % Assign pattern name
            try
                Patterns(iPartition,j,i).name = Option.patternNames(i);
            catch 
                Patterns(iPartition,j,i).name = "Pattern " + i;
            end

            if iPartition == 1
                if source == "hpc"
                    s_all = r.hpc.X{i};
                    s_ind_all = 1:size(s_all, 1);
                elseif source == "pfc"
                    s_all = r.pfc.X{i};
                    s_ind_all = 1:size(s_all, 1);
                else
                    error('Source area not recognized')
                end
                if target == "hpc"
                    t_all = r.hpc.X{i};
                    t_ind_all = 1:size(t_all, 1);
                elseif target == "pfc"
                    t_all = r.pfc.X{i};
                    t_ind_all = 1:size(t_all, 1);
                else
                    error('Target area not recognized')
                end
                if source ~= target
                    assert(numel(s_ind_all) ~= numel(t_ind_all), 'Source and target areas are the same')
                end
                % Assign x_source and x_target
                Patterns_overall(j,i).X_source = s_all;
                Patterns_overall(j,i).X_target = t_all;
                % Assign index_source and index_target
                Patterns_overall(j,i).index_source = s_ind_all;
                Patterns_overall(j,i).index_target = t_ind_all;
                % Assign directionality
                Patterns_overall(j,i).directionality = directionality;
                % Assign pattern name
                try
                    Patterns_overall(j,i).name = Option.patternNames(i);
                catch 
                    Patterns_overall(j,i).name = "Pattern " + i;
                end
            end

        end
    end
end

% Add overall for ALL TIME
areaPerNeuron = r.areaPerNeuron;
sessionTypePerBin = r.sessionTypePerBin;
for j = 1:numel(split_info.directionality)
    if i == HPC
        s_all = r.spikeRateMatrix(:, sessionTypePerBin == 1 & ... 
            areaPerNeuron == "CA1");
        t_all = s_all;
        s_inds = r.celllookup(areaPerNeuron == "CA1").index_by_region;
        t_inds = s_inds;
        source = "hpc";
        target = "hpc";
    else
        s_all = r.spikeRateMatrix(:, sessionTypePerBin == 1 & ... 
            areaPerNeuron == "CA1");
        t_all = r.spikeRateMatrix(:, sessionTypePerBin == 1 & ... 
            areaPerNeuron == "PFC")
        s_inds = r.celllookup(areaPerNeuron == "CA1").index_by_region;
        t_inds = r.celllookup(areaPerNeuron == "PFC").index_by_region;
        source = "hpc";
        target = "pfc";
    end
    % Assign x_source and x_target
    Patterns_overall(j, end+1).X_source = s_all;
    Patterns_overall(j, end).X_target = t_all;
    % Assign index_source and index_target
    Patterns_overall(j, end).index_source = s_inds;
    Patterns_overall(j, end).index_target = t_inds;
    % Assign directionality
    Patterns_overall(j, end).directionality = directionality;
    % Assign pattern name
    Patterns_overall(j, end).name = "Overall";
    time = r.timeBinMidPoints(sessionTypePerBin == 1);
    Patterns_overall(j, end).X_time = reshape(time',1,[]);
    Patterns_overall(j, end).source = source;
    Patterns_overall(j, end).target = target;
end

szCellOfWindows = squeeze(size(r.windowInfo.cellOfWindows));
if numel(szCellOfWindows) == 2 && szCellOfWindows(1) == 1
    szCellOfWindows = szCellOfWindows(2);
end

Patterns = reshape(Patterns, ...
    [Option.numPartition, numel(split_info.directionality), ...
    size(r.windowInfo.cellOfWindows)]);

disp(['Partitioning data took ', num2str(toc), ' seconds.'])


