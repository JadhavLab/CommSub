%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   
% ,---.o         |    |              |        |     
% |__. .,---.,---|    |    ,---.,-.-.|---.,---|,---.
% |    ||   ||   |    |    ,---|| | ||   ||   |,---|
% `    ``   '`---'    `---'`---^` ' '`---'`---'`---^
%                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script finds a good lambda for the spectral 
% data.

% Parameters
% ----------
bands               = [6 12; 0.5 4]; % 1- theta, 2 - delta
lowhigh_defintions = [0.2, 0.8]; % percentile definitions for high and low respectively
%animal_list         = {'JS12','JS13','JS14','JS15','JS21','ER1','KL8'};
animal_list         = {'JS15','JS21','ER1','KL8'};
fields              = {'S1','S2','wpli'};
CA1 = 1; PFC = 2; CA1PFC = 3;
theta = 1; delta = 2;
high = 2; low = 1;
frequency = 2; time = 1;
samprate = 0.1;
win = [10, 10]; % 10 seconds before and after
%plot_range = [0, 40];



% Add animal data to path
% ------------------------
for animal = animal_list
    animalToPath(animal{1});
end

% figures out optimal lambda over the dataset

% Tracked
%lambdas = sort([logspace(-1,-5,nLambdas)], 'ascend'); 
%nLambdas = 100; % increase if you're patient
%lambdas = sort([logspace(-3,-5,nLambdas)], 'ascend'); 
%loadings       = [];
%regularization = [];
%cost           = [];
%lambda_optimi  = [];

% Iterate each animal and run seqnmf logic
% ----------------------------------------
for animal = animal_list


    % Determine the days and epochs for the animal
    animaldef_ = animaldef(animal{1});
    task = loaddatastruct(animaldef_{2:3}, 'task');
    search_results = cellfetch(task, 'environment');
    dayepochs = search_results.index;


    % Load up all data for an animal
    % ------------------------------
    P = []; % for aggreagting pos
    C = []; % for aggreagting pos
    for dayepoch = dayepochs'

        fprintf('Day %d, epoch %d of animal %s\n', dayepoch, animal{1})

        % -- load spectral info --
        cgramc_file = [animal{1} 'cgramc' sprintf('-%02d-%02d.mat', dayepoch)];
        if ~exist(cgramc_file, 'file')
            continue
        end
        load(cgramc_file)
        data = ffend(cgramc);
        % -- load position info --
        load([animal{1} 'pos' sprintf('%02d', dayepoch(1))])
        day   = dayepoch(1);
        epoch = dayepoch(2);
        p = pos{day}{epoch};
        cgramc = ffend(cgramc);
        try
            C = [C; cgramc];
        catch ME
            warning(sprintf('Day %d, epoch %d of animal %s FAILED\n', dayepoch, animal{1}))
        end
    end

    % Aggregate
    % ---------
    data.t = cat(2, C.t);
    data.f = cat(1, C.f);
    for field = fields
        data.(field{1}) = cat(1, C.(field{1}));
    end



    data.data = seqnmf_packfields(data, fields, 'zscore', fields(1:2), 'minmax', fields(3), 'floorceil', [-1,2])';
    X = data.data;
    X(isnan(X)) = 0;
    X = X(:,end-round(5*mean(arrayfun(@(x) numel(x.t), C))):end); % use 5 epochs of data

    %,---.|                   o              |              |        |     
    %|    |---.,---.,---.,---..,---.,---.    |    ,---.,-.-.|---.,---|,---.
    %|    |   ||   ||   |`---.||   ||   |    |    ,---|| | ||   ||   |,---|
    %`---'`   '`---'`---'`---'``   '`---|    `---'`---^` ' '`---'`---'`---^
    %                               `---'                                  
    %% Procedure for choosing lambda
    K = 10; 
    Lneural = round(10/cgramc.params.movingwin(2)); % Replace this with optomized L value
    for li = 1:length(lambdas)
        [N,T] = size(X);
        tic
        % 'lambdaL1W', .1,
        [W, H, ~, loading, power]= seqNMF_gpu(X,'K',K,'L',Lneural,...
             'lambda', lambdas(li), 'maxiter', 50, 'showPlot', 0); 
        toc
        assert(~all(isnan(W(:))));
        loadings(end+1,:) = loading;
        [cost(end+1), regularization(end+1),~] = helper.get_seqNMF_cost(X,W,H);
        display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
        LD = repmat(lambdas, [1, ceil(numel(cost)/numel(lambdas))]);
        sfigure(50); clf; p1=plot(LD(1:numel(cost)), cost, 'o-'); hold on; p2=plot(LD(1:numel(cost)), regularization, 'o-'); set(gca,'yscale','log'); legend([p1,p2],'Cost','Regularization'); ylim([0,inf]);
        numel(regularization)
    end

    % checkpointing
    if ~exist('~/Data/deltacoherence/findlambda/', 'dir')
        mkdir('~/Data/deltacoherence/findlambda/')
    end
    save('~/Data/deltacoherence/findlambda/checkpoint.mat', 'lambdas', 'loadings', 'regularization', 'cost');

end

LD = repmat(lambdas, [1, numel(cost)/numel(lambdas)]);
animals = []
for animal = animal_list
    animals = [animals, repmat(string(animal{1}), [1, numel(lambdas)])];
end

[~,idx] = sort(LD);
LD = LD(idx);
R = regularization(idx);
C = cost(idx);
C = reshape(C', 3, [])';

%make lambdas equal in length to the number of cost iterations
%% plot costs as a function of lambda
windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
Rs = filtfilt(b,a,R); 
minRs = prctile(R,10); maxRs= prctile(R,90);
Rs = (Rs-minRs)  /  (maxRs-minRs); 
R =  (R-minRs)   /  (maxRs-minRs); 
minCs =  prctile(C,10); maxCs =  prctile(C,90); 
Cs = filtfilt(b,a,C); 
C = (C -minCs)./(maxCs-minCs); 
c = C';
C = c(:);
Cs = (Cs -minCs)./(maxCs-minCs); 
c = Cs';
Cs = c(:);

figure(51);
clf; 
plot(LD,Rs, 'b');hold on
plot(LD,Cs,'r')
scatter(LD, R, 'b', 'markerfacecolor', 'flat');
scatter(LD, C, 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost (au)')
set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
%
%
