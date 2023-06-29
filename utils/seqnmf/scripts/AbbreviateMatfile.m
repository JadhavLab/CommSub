allcandidates = false;

% The name of the root folder
rootFolder = '/Volumes/Data/deltacoherence/';

fieldsToPreserve = {'t',  'f', 'param', 'date', 'animalnames', 'animalcnt', 'epoch', 'inds', 'subsample', 'seqStyle', 'H', 'K', 'L', 'WH', 'W', ...
                    'cost','loadings','power','samprate','fields','subsample_indices', 'subsample_logical', 'behavior2data_inds', 'data'};

% Replots all of the results from SeqNMF
InitializeParamsets;
% Iterates over each of the parameter sets and replots
% the data
diary(['~/Projects/deltacoherence/results/diary__' date])
for param = progress( paramsets, 'Title', 'Parameter iteration')

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
    data = M.data;
    fieldsToRemove = [];
    numericFields = [];
    for field = fieldnames(data)'
        field = field{1};
        if ~(ismember(field, fieldsToPreserve) || contains(field, data.fields))
            fieldsToRemove = [fieldsToRemove, string(field)];
        elseif isnumeric(data.(field))
            data.(field) = single(data.(field));
        end
    end
    fieldsToRemove
    if all(ismember({'WHprime','Wprime'}, fieldnames(data)))
        M.Properties.Writable = true;
        M.data = rmfield(data,{'WHprime','Wprime'});
    end
    behavior = M.behavior;

    % Create exclusive pre file with data 
    % Make abbreviated matfile
    clear M;
    data = rmfield(data, fieldsToRemove);
    disp("Wrote file with the following fields:")
    disp(string(sprintf('\n')) + fieldnames(data));
    M.behavior = behavior;
end
