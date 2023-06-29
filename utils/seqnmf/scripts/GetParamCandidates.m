% The name of the root folder
rootFolder = '/Volumes/Data/deltacoherence/';

% Replots all of the results from SeqNMF
InitializeParamsets;
% Iterates over each of the parameter sets and replots
% the data
diary(['~/Projects/deltacoherence/results/diary__' date])
ParamCandidates = {};
for param = progress( paramsets(1:end), 'Title', 'Parameter iteration')

    % Deal out each of the parameters into unique variables
    [epoch_type, timescale, K, fieldstr, seqStyle] = deal( param{1}{:} );
    % Find the most recent matching folder
    % ------------------------------------
    seq.initialize('epoch_type', epoch_type,'timescale', timescale, 'fieldstr',...
                    fieldstr, 'orthoName', seqStyle, 'skipList', {'K'}, 'K', K);
    %eqnmf_kws{1}{end}=7; seqnmf_kws{1}{2} = 1e-4; seqnmf_kws{1}{6} = 4e-4; seqnmf_kws{1}{4} = 4e-4;
    kws = seq.container2kws(seqnmf_kws);
    GetRecentMatfile;
    if notFoundFlag || missingFlag
        continue
    end
end
