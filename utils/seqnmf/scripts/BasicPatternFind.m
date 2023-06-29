% ===========================================================================
% -------------
% Things to do
% -------------
% 1. optomize the lambda
% 2. optomize L? (lower priority, seems to remove patterns on its own --
% possibly speeds)
% 3. Another animal
%   a. Correlate
%   b. W(animal_1) = W(animal_2) in terms of loading
% 4. Spectral frequencies of a W?
% ===========================================================================

% Flags
skipMaster = true;
storeDatBeforeSubsamp = true;
subsample = 0.09; % fraction of data to sample

% -----------------------
% Add animal data to path
% ------------------------
animal_list = {'JS12','JS13','JS14','JS15','JS21','ER1','KL8'};
for animal = animal_list
    animalToPath(animal{1});
end

InitializeParamsets;

%% Preamble
diary(['~/Projects/deltacoherence/results/diary_' date])
for param = progress( paramsets(1:end), 'Title', 'Parameter iteration')
    

    % Deal out each of the parameters into unique variables
    [epoch_type, timescale, K, fieldstr, seqStyle] = deal( param{1}{:} );

    % ----------------------------------------------------------
    %     _                   _                _       _        
    %    / \   ___ __ _ _   _(_)_ __ ___    __| | __ _| |_ __ _ 
    %   / _ \ / __/ _` | | | | | '__/ _ \  / _` |/ _` | __/ _` |
    %  / ___ \ (_| (_| | |_| | | | |  __/ | (_| | (_| | || (_| |
    % /_/   \_\___\__, |\__,_|_|_|  \___|  \__,_|\__,_|\__\__,_|
    %                |_|                                        
    % ----------------------------------------------------------

    % Initialize folders and parameters
    % ---------------------------------
    seq.initialize('K', K, 'epoch_type', epoch_type,'timescale', timescale, 'fieldstr',...
                    fieldstr, 'orthoName', seqStyle, 'skipList', {'K'}, 'maxiter', 450);
    %eqnmf_kws{1}{end}=7; seqnmf_kws{1}{2} = 1e-4; seqnmf_kws{1}{6} = 4e-4; seqnmf_kws{1}{4} = 4e-4;
    kws = seq.container2kws(seqnmf_kws);
    folder = seqnmf_folder(fieldstr, kws{:}, 'timescale', timescale, 'epoch_type', char(epoch_type));

    % Announce the folder and create it
    % ---------------------------------
    disp("Starting " + string(seqStyle) + " in folder " + string(folder))
    disp(sprintf('epoch=%s timescale=%d K=%d fieldstr=%s seqstyle=%s', param{1}{:}))
    assert(~(ismember('animal_overall', nmfs) && ismember('overall', nmfs)))
    fullfolder = sprintf('~/Data/deltacoherence/%s/%s', folder, seqStyle);
    if ~exist(sprintf('~/Data/deltacoherence/%s', folder), 'dir')
        mkdir(sprintf('~/Data/deltacoherence/%s', folder))
    end
    if ~exist(fullfolder, 'dir')
        mkdir(fullfolder)
    end
    file = fullfile(fullfolder, 'master_seqnmf')

    GetData; % Relies on animal_list existing
    
    RunSeqnmf; % Run seqnmf

    %% MASTER DATASET
    % Efizz dataset
    M = matfile(file, 'Writable', true);
    if ~skipMaster
        D = {data.t', data.data'};
        D = table(D{:}, 'VariableNames', {'time','data'});
        % Master dataset
        if isequal(epoch_type, 'run')
            inds = lookup(double(D.time), double(behavior.time));
            master = [behavior(inds,:), table(D.data,'VariableNames',{'data'})];
        else
            master = D;
        end
        % Eliminate large time skips
        master.time(:) = [0; diff(master.time(:))];
        master.time(master.time(:) > 10) = 0;
        master.time(:) = cumsum(master.time(:));
        M.master = master;
        clear D
    end
    %% MASTER DATASET LIGHT
    % Efizz dataset
    %M = matfile(replace(file,'nmf','nmf_light'), 'Writable', true);
    M = matfile(replace(file,'nmf','nmf'), 'Writable', true);
    M.behavior = behavior;
    M.data = rmfield(data,'pre');
    M.pre = data.pre;
    clear M;
    %M.behavior = behavior;
    %M.data = data;
    %clear M;

    % -----------------------------
    %  ____  _     ___ _____ ____  
    % |  _ \| |   / _ \_   _/ ___| 
    % | |_) | |  | | | || | \___ \ 
    % |  __/| |__| |_| || |  ___) |
    % |_|   |_____\___/ |_| |____/ 
    % -----------------------------

    if ~skipMaster
        % RAW DATA
        data_plot = fig('rawData'); clf
        fstring = string(F);
        file = sprintf('~/Data/deltacoherence/%s/data', epoch_type);
        mkdir(fileparts(file))
        characterize.rawData(master, epoch_type);
        % Title
        suptitle(sprintf('Data: ArgNum: %s, Timescale: %d, fields: %s', seqStyle, timescale, fstring.join('-')));
        saveas(data_plot, sprintf('%s.%s', file, 'svg'))
        try
            saveas(data_plot, sprintf('%s.%s', file, 'fig'))
        catch ME
            save(sprintf('%s.%s', file, 'fig.mat'), 'data_plot', '-v7.3')
        end
    end

    % Loadings
    loading_plot = fig(char(fullfile('~/Projects/deltacoherence/results', folder, seqStyle, 'loadings')));
    plot(1:K, data.loadings);
    saveas(loading_plot, [loading_plot.Name '.svg'])
    saveas(loading_plot, [loading_plot.Name '.fig'])

end
%end % seqStyle
%end % timescale
%end % fields
%end % epoch_type (run, sleep)
