function Patterns = partition(r, Option) 
% ----------------------------
% Note:
% ----------------------------
% when the partition is three-ways, direction==1 means same target/source pair
% and direction==2 means diff target/source pair
% ----------------------------
 
% Initialize the scaffold of the pattern struct
Patterns = initPatternStruct();

% Split cells into source and target areas
[s_hpc, s_pfc, t_hpc, t_pfc, s_hpc_index,...
 s_pfc_index, t_hpc_index, t_pfc_index, directionality] = ...
    trialSpikes.split.legacyRun(r, Option);

% ------------------------------
% Place paritioned data properly
% ------------------------------
for iPartition = 1:Option.numPartition
for i = 1:numel(r.windowInfo.cellOfWindows)
    for j = 1:numel(directionality)
        
        % Parse directionality
        sourcetarg = directionality.split('-');
        source = sourcetarg(1);
        target = sourcetarg(1);
        if source == "hpc"
            s_dat = s_hpc;
            s_ind = s_hpc_index;
        elseif source == "pfc"
            s_dat = s_pfc;
            s_ind = s_pfc_index;
        end
        if target == "hpc"
            t_dat = t_hpc;
            t_ind = t_hpc_index;
        elseif target == "pfc"
            t_dat = t_pfc;
            t_ind = t_pfc_index;
        end
        
        % Assign x_source and x_target
        Patterns(iPartition,j,i).X_source = s_dat{i};
        Patterns(iPartition,j,i).X_target = t_dat{i};
        
        % Assign index_source and index_target
        Patterns(iPartition,j,i).index_source = s_ind;
        Patterns(iPartition,j,i).index_target = t_ind;
        
        % Assign directionality
        Patterns(iPartition,j,i).directionality = directionality(j);
        
        % Assign pattern name
        try
            Patterns(iPartition,j,i).name = Option.patternNames(i);
        catch E
        end
    end
end
end

szCellOfWindows = squeeze(size(r.windowInfo.cellOfWindows));
if numel(szCellOfWindows) == 2 && szCellOfWindows(1) == 1
    szCellOfWindows = szCellOfWindows(2);
end
Patterns = reshape(Patterns, [Option.numPartition, numel(directionality), size(r.windowInfo.cellOfWindows)]);

