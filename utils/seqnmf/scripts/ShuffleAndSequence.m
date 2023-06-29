% ------------------------
% Boilerplate
% --------------------
seqnmf_initialize; % Boilerplate across analyses
% Lower number of iterations to minimal amount for good convergence
% Set pattern to number found and remove regularization that is used to wittle that pattern downw.
for i = 1:numel(seqnmf_kws)
    f = find(cellfun(@(x) isequal(x, 'maxiter'), seqnmf_kws{i}));
    seqnmf_kws{i}{f+1} = 100; % According to cost function progression, more than enough iterations to asymptotically minimize
    %f = find(cellfun(@(x) isequal(x, 'lambda'), seqnmf_kws{i}));
    %seqnmf_kws{i}{f+1} = 0; % For checking shuffles and sequenciness, remmove the regularization pressure
    f = find(cellfun(@(x) isequal(x, 'K'), seqnmf_kws{i}));
    seqnmf_kws{i}{f+1} = 3; % For my purposes, bring K to the actual value found every time I run this, which is 1
end
folder = [folder filesep 'sequenciness+shuffle'];
disp(['Folder = ' folder]);

% How many iterations to perform
% ------------------------------
nIter         = 1;   % How many to maximize over per shuffle
nRepsShuff    = 100; % Number of raw shuffles
nRepsColShuff = 100; % Number of column shuffles

% ------------------ TIME ESTIMATE -----------------------------
% TIME ESIMATE: 16.81 seconds per iteration with a serious GPU
% (30*(100*5 + 100*5) * one_iter)/3600)/24 = 5.8 Days
% --------------------------------------------------------------

% ------------------------
% Add animal data to path
% ------------------------
for animal = animal_list
    animalToPath(animal{1});
end

% ----------------------------------------
% Iterate each animal and run seqnmf logic
% ----------------------------------------
disp("Starting " + string(seqnmf_arg_names(argNum)))
P = []; % for aggreagting pos
C = []; % for aggreagting pos
A = []
for animal = animal_list


    % Determine the days and epochs for the animal
    animaldef_ = animaldef(animal{1});
    task = loaddatastruct(animaldef_{2:3}, 'task');
    search_results = cellfetch(task, 'environment');
    dayepochs = search_results.index;


    % Mode-specific behaviors
    if ismember('animal_overall', nmfs)
        disp('Clearing for next animal')
        P = []; % for aggreagting pos
        C = []; % for aggreagting pos
        file = @(x) sprintf('~/Data/deltacoherence/%s/%s/%sseqnmf_repi=%d.mat', folder, seqnmf_arg_names(argNum), animal{1}, x);
        if exist(file,'file')
            disp('Continuing')
            continue
            fprintf('Redoing %s...\n', file);
        else
            fprintf('Working on %s...\n', file);
        end
    end

    % Load up all data for an animal
    % ------------------------------
    for dayepoch = dayepochs'
        % -- load spectral info --
        cgramc_file = [animal{1} 'cgramc' sprintf('-%02d-%02d.mat', dayepoch)];
        if ~exist(cgramc_file, 'file')
            disp('Continuing')
            continue
        else
            disp("Processing")
        end
        load(cgramc_file)
        data = ffend(cgramc);
        % -- load position info --
        load([animal{1} 'pos' sprintf('%02d', dayepoch(1))])
        day   = dayepoch(1);
        epoch = dayepoch(2);
        p = pos{day}{epoch};
        P = [P; p.data];
        try
            C = [C; data];
        catch ME
            warning(sprintf('Skipping animal %s, day %d epoch %d\n', animal{1}, dayepoch));
        end
        L =  round(timescale/data.params.movingwin(2));

        if ismember('epoch', nmfs)
            data  = ry_run_seqnmf(data, fields, 'seqnmf_kws', {seqnmf_kws{argNum}{:}, 'L', L, 'checkpoint', 200}, 'fieldpack_kws', fieldpack_kws);
        end
    end
end

% -----------
% OVERALL
% -----------
if ismember('overall', nmfs)
        file = @(x) sprintf('~/Data/deltacoherence/%s/%s/master_seqnmf_rep=%d.mat', folder, seqnmf_arg_names(argNum), x);
        % Aggregate
        % ---------
        data.animal = cell(1,numel({C.t})); 
        cnt = 0;
        for t = {C.t}
            cnt = cnt + 1;
            data.animal{cnt} = cnt * ones(size(t{1}));
        end
        data.t = cat(2, C.t);
        data.animal = cat(2, data.animal{:});
        for field = fields
            data.(field{1}) = cat(1, C.(field{1}));
        end
        F = fields;
        if any(contains(fields, 'phi'))
            data.phi_sin = sin(data.phi);
            data.phi_cos = cos(data.phi);
            F{cellfun(@(x) isequal(x, 'phi'), F)} = 'phi_sin';
            F = [F 'phi_cos'];
        end

        fieldpack_tmp = [fieldpack_kws, 'groups', data.animal];
        data.data = seqnmf_packfields(data, fields, fieldpack_tmp{:}) ;
        X = data.data';
        PEx = [];
        PExColShuff = [];
        PExShuff = [];

        % ---------------------
        tic
        disp('----------------------')
        disp('--- SIGNAL -----------')
        disp('----------------------')
        for iteri = 1:1
            disp(['Iteration ' num2str(iteri)])
            rng('shuffle')
            seqnmf_factorize = @seqNMF_gpu;
            [W, H, ~,~,tmp] = seqnmf_factorize(X, seqnmf_kws{argNum}{:}, 'showPlot', 0);
            PEx(iteri) = tmp;
            % Save one example
            if iteri == 2
                save(sprintf('~/Data/deltacoherence/%s/PEx_example.mat',folder), 'W', 'H', 'X')
            end
        end
        disp('Signal finished')
        toc
        save(sprintf('~/Data/deltacoherence/%s/shuffleseq.mat',folder), 'PEx')

        % Run on TOTAL shuffle
        % ---------------------
        [N T] = size(X)
        disp('----------------------')
        disp('--- SHUFFLE ----------')
        disp('----------------------')
        for repi = 88:nRepsShuff
            tic
            disp(['Replication ' num2str(repi)])
            Xshuff = X;
            for ni = 1:N
                timeshuff = randperm(T);
                Xshuff(ni,:) = X(ni, timeshuff); 
            end
            
            for iteri = 1:nIter
                rng('shuffle')
                seqnmf_factorize = @seqNMF_gpu;
                [W, H, ~,~,PExShuff(iteri,repi)] = seqNMF_gpu(Xshuff, seqnmf_kws{argNum}{:}, 'showPlot',0);
            end
            % Save one example
            if repi == 1
                save(sprintf('~/Data/deltacoherence/%s/PExShuff_example.mat',folder), 'W', 'H', 'X', 'Xshuff')
            end
            toc
        end
        save(sprintf('~/Data/deltacoherence/%s/shuffleseq.mat',folder), 'PEx', 'PExShuff')
        % Left off on repi=22

        % Run on COLUMN shuffle
        % ---------------------
        for repi = 9:nRepsColShuff
            disp(['Replication ' num2str(repi)])
            Xshuff = X(:,[1:L (L + randperm(T-L))]); % don't shuffle to first L bins... these cannot be explained by seqNMF
            tmp = [];
            
            for iteri = 1:nIter
                rng('shuffle')
                [W, H, ~,~,PExColShuff(iteri, repi)] = seqNMF_gpu(X, seqnmf_kws{argNum}{:}, 'showPlot', 0);
            end    
            % Save one example
            if repi == 1
                save(sprintf('~/Data/deltacoherence/%s/PExColShuff.mat',folder), 'W', 'H', 'X', 'Xshuff')
            end
        end
        % Repi = 22
        save(sprintf('~/Data/deltacoherence/%s/shuffleseq.mat',folder), 'PEx', 'PExShuff','PExColShuff')

        % Obtain sequenciness statistic
        % -----------------------------
        NoiseFloor  = median(max(PExShuff,[],1));
        NoiseFloor  = (max(PExShuff,[],1));
        SyncFloor   = median(max(PExColShuff,[],1));
        SyncFloor   = (max(PExColShuff,[],1));
        TotalSignal = max(PEx(:),[], 1);
        PAS         = (TotalSignal - SyncFloor)./...
                      (TotalSignal - NoiseFloor)
        
        save(sprintf('~/Projects/deltacoherence/results/%s/%s/master_shuffle.svg', folder, seqnmf_arg_names(argNum)), 'PEx','PExShuff','PExColShuff','NoiseFloor','SyncFloor','TotalSignal')

        sfigure;
        subplot(1,2,1)
        b=[];
        b(1) = bar(1, TotalSignal);
        b(2) = bar(2, SyncFloor);
        b(3) = bar(3, NoiseFloor);
        text(0.25, 0.6, 'Sequenciness score = %2.2f', PAS);
        legend(b, 'Total Signal', 'Synchronous Signal Floor', 'Noise Floor')
        subplot(1,2,2)
        h=[];
        h(1) = histogram(TotalSignal(:));
        h(2) = histogram(SyncFloor(:));
        h(3) = histogram(NoiseFloor(:));
        legend(h, 'Total Signal', 'Synchronous Signal Floor', 'Noise Floor')

        saveas(gcf, sprintf('~/Projects/deltacoherence/results/%s/%s/master_shuffle.svg', folder, seqnmf_arg_names(argNum)));
        saveas(gcf, sprintf('~/Projects/deltacoherence/results/%s/%s/master_shuffle.fig', folder, seqnmf_arg_names(argNum)));

end
