% Flags
skipMaster = true;
storeDatBeforeSubsamp = true;
subsample = 0.085; % fraction of data to sample
resultsFolder = project('seqNMF');

% -----------------------
% Add animal data to path
% ------------------------
animal_list = {'JS12','JS13','JS14','JS15','JS21','ER1','KL8'};
for animal = animal_list
    animalToPath(animal{1});
end

InitializeParamsets;
disp('Available paramsets:')
paramsetstable = unique(removevars(paramsetstable,'K'))

% Deal out each of the parameters into unique variables
disp('Picked set:')
pick = paramsetstable(1,:),
params = table2cell(pick);
[epoch_type, timescale, fieldstr, seqStyle] = deal( params{:} );
seq.initialize('epoch_type', epoch_type,'timescale', timescale, 'fieldstr',...
                fieldstr, 'orthoName', seqStyle, 'skipList', {'K'}, 'maxiter', 450);
seqnmf_kws.remove({'K','lambda','lambdaOrthoH','lambdaOrthoW'})
kws = seq.container2kws(seqnmf_kws);
folder = seqnmf_folder(fieldstr, kws{:}, 'timescale', timescale, ...
                        'epoch_type', char(epoch_type), 'usedate', false);
GetData;

X = data.data;
trainNEURAL = data.data;
testNEURAL = data.pre.data(:,~data.subsample_logical);

mkdir(fullfile(resultsFolder,'findParams'));
datfolder=fullfile(resultsFolder,'findParams', folder);
if exist(datfolder,'dir')
    disp('Result exists');
    keyboard;
end
mfile = matfile(fullfile(resultsFolder, 'findParams', 'choosek' + string(folder)), 'Writable', true);
mfile.params = pick;

%
%                                              
% ,---.|                   o              |   /
% |    |---.,---.,---.,---..,---.,---.    |__/ 
% |    |   ||   ||   |`---.||   ||   |    |  \ 
% `---'`   '`---'`---'`---'``   '`---|    `   `
%                                `---'         
%% Procedure for choosing K
tic
mfile.data = data;
Ws = {};
Hs = {};
numfits = 3; %number of fits to compare
kset = 1:12;
ksubset = kset;

p = ProgressBar(numel(ksubset), 'Title', 'Choosing K')
for k = ksubset
    if false
    display(sprintf('running seqNMF with K = %i',k))
    parfor ii = 1:numfits
        [Ws{ii,k},Hs{ii,k}] = seqNMF(X,'K',k, 'L', L,'lambda', 0,'maxiter',140,'showplot',0); 
        % note that max iter set low (30iter) for speed in demo (not recommended in practice)
    end

    end
    inds = nchoosek(1:numfits,2);
    for i = 1:size(inds,1) % consider using parfor for larger numfits
            Diss(i,k) = helper.DISSX(Hs{inds(i,1),k},Ws{inds(i,1),k},Hs{inds(i,2),k},Ws{inds(i,2),k});
    end
    choosek = struct('Diss',Diss,'Ws',Ws,'Hs',Hs);
    mfile.choosek = choosek;
    p.step([],[],[]);
end
toc

%% Plot Diss and choose K with the minimum average diss.
f=fig("Choose K - " + folder)
plot(kset, Diss,'ko'), hold on
h1 = plot(kset, median(Diss,1),'k-','linewidth',2);
[~,K]  = max(diff(median(Diss,1)))
h2 = plot([K,K],[0,0.5],'r--');
legend([h1 h2], {'Median Diss','True K'})
xlabel('K')
ylabel('Diss')
choosek.f = f;
choosek.K = K;
mfile.choosek = choosek


% display('Testing significance of factors on held-out data')
% [pvals,is_significant] = test_significance(testNEURAL,W)
% 
% W = W(:,is_significant,:); 
% H = H(is_significant,:); 

%% plot, sorting neurons by latency within each factor
%[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
%indSort = hybrid(:,3);
%tstart = 180; % plot data starting at this timebin
%figure; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), ...
%    0,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end))
%title('Significant seqNMF factors, with raw data')
%figure; WHPlot(W(indSort,:,:),H(:,tstart:end), ...
%    helper.reconstruct(W(indSort,:,:),H(:,tstart:end)),...
%    0,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end))
%title('SeqNMF reconstruction')


% ,---.|                   o              |              |        |     
% |    |---.,---.,---.,---..,---.,---.    |    ,---.,-.-.|---.,---|,---.
% |    |   ||   ||   |`---.||   ||   |    |    ,---|| | ||   ||   |,---|
% `---'`   '`---'`---'`---'``   '`---|    `---'`---^` ' '`---'`---'`---^
%                                `---'                                  
%% Procedure for choosing lambda
nLambdas = 30; % increase if you're patient
X = trainNEURAL;
lambdas = sort([logspace(-1,-6,nLambdas)], 'ascend'); 
W1penalties = [0, 1e-3, 1e-2,1e-1];
loadings = [];
cost  = zeros(length(W1penalties),length(lambdas));
regularization  = zeros(length(W1penalties),length(lambdas));
[Hs, Ws] = deal({});
power  = zeros(length(W1penalties),length(lambdas));
clear chooseLambda
chooseLambda = struct();
chooseLambda.Ws = cell();
chooseLambda.Ws = cell();

prog = ProgressBar(numel(W1penalties) * numel(lambdas), 'Title', 'Choosing lambda')
for li = 1:numel(lambdas)
    for wi = 1:numel(W1penalties)
        W1penality = W1penalties(wi);
        [N,T] = size(X);
        %try
            [Ws{wi,li}, Hs{wi,li}, C,loadings(wi, li,:),power(wi,li)] = ...
                seqNMF_gpu(X,'K',single(K),'L',single(L),...
                'lambda', single(lambdas(li)), 'maxiter', 150, 'showPlot', 0, 'lambdaL1H', single(W1penality)); 
            %[Ws{wi,li}, Hs{wi,li}, C,loadings(wi, li,:),power(wi,li)] = ...
            %    seqNMF(X,'K',single(K),'L',single(L),...
            %    'lambda', single(lambdas(li)), 'maxiter', 150, 'showPlot', 0, 'lambdaL1H', single(W1penality)); 
            prog.step([],[],[]);
            [cost(wi, li),regularization(wi, li),~] = helper.get_seqNMF_cost(X,Ws{wi,li},Hs{wi,li});
            display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])

            %li = li + 1;
        %catch ME
        %    disp('Failed, restarting')
        %    keyboard
        %    reset(gpuDevice)
        %end
    end

    vals = {'loadings', loadings, 'power', power, ...
                            'cost', cost, 'regularization', regularization, ...
                            'W1penalties', W1penalties, 'lambdas', lambdas};
    for v = 1:2:numel(vals)
        chooseLambda.(vals{v}) = vals{v+1};
    end
    chooseLambda.Hs = Hs;
    chooseLambda.Ws = Ws;

    mfile.chooseLambda = chooseLambda;
end

lambdainterp = linspace(min(lambdas), max(lambdas), 5e6);
I = @(x,y) interp1(x,y,lambdainterp,'linear');

clear F
F= {};
targetLambdas = [];
%% plot costs as a function of lambda
for wi = 1:length(W1penalties)
    fig([folder ' WL1 = ' + string(W1penalties(wi))]); clf
    F{end+1} = gcf;
    windowSize = 3; 
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    %Filter regularization penalty
    Rs = filtfilt(b,a,regularization(wi,:)); % Filtered regularization
    % Normalize reg penalty
    minRs = prctile(regularization(wi,:),10); 
    maxRs = prctile(regularization(wi,:),90);
    Rs = (Rs-minRs)/(maxRs-minRs); % filtered: 0 to 1
    R = (regularization(wi,:)-minRs)/(maxRs-minRs);  % unfiltered: 0 to 1
    Cs = filtfilt(b,a,cost(wi,:));  % filter cost
    minCs =  prctile(cost(wi,:),10); maxCs =  prctile(cost(wi,:),90);
    Cs = (Cs -minCs)/(maxCs-minCs); % filtered : Scale cost 0 to 1
    C = (cost(wi,:) -minCs)/(maxCs-minCs);  % unfiltered : Scale cost 0 to 1
    % Plot crossing
    clf; hold on
    plot(lambdainterp,I(lambdas,Rs), 'b')
    plot(lambdainterp,I(lambdas,Cs),'r')
    scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
    scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
    xlabel('Lambda'); ylabel('Cost (au)')
    set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
    set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
    set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
    [~,targetLambdas(end+1)] = min( abs(I(lambdas, Cs)-I(lambdas,Rs)) );
    targetLambdas(end) = lambdainterp(targetLambdas(end));
end

fig([folder ' W1 relation to penalty']); clf;
F{end+1} = gcf;
plot(W1penalties, targetLambdas, '*:');
set(gca,'xscale','linear');
title(sprintf('Corrcoeff = %0.2f %0.2f %0.2f %0.2f', corrcoef([W1penalties', targetLambdas'])));
chooseLambda = struct('loadings', loadings, 'power', power, ...
                        'cost', cost, 'regularization', regularization, ...
                        'W1penalties', W1penalties, 'lambdas', lambdas, 'targetLambdas', targetLambdas);
chooseLambda.Ws = Ws;
chooseLambda.Hs = Hs;
chooseLambda.F = F;
mfile.chooseLambda = chooseLambda;


%Code to fix improper save structure
chooseLambda = mfile.chooseLambda;
Ws = reshape({chooseLambda.Ws}, size(chooseLambda));
Hs = reshape({chooseLambda.Hs}, size(chooseLambda));
chooseLambda = chooseLambda(end);
chooseLambda.W1penalties = W1penalties;
chooseLambda.lambdas     = targetLambdas;
chooseLambda.Ws     = Ws;
chooseLambda.Hs     = Hs;
mfile.chooseLambda = chooseLambda;


%                                    
%,---.|                   o          
%|    |---.,---.,---.,---..,---.,---.
%|    |   ||   ||   |`---.||   ||   |
%`---'`   '`---'`---'`---'``   '`---|
%                               `---'
%                                                            
%          |    |         |   ||              |        |     
%,---.,---.|--- |---.,---.|---||    ,---.,-.-.|---.,---|,---.
%|   ||    |    |   ||   ||   ||    ,---|| | ||   ||   |,---|
%`---'`    `---'`   '`---'`   '`---'`---^` ' '`---'`---'`---^
chooseLambda = mfile.chooseLambda;

%% Procedure for choosing lambda
paramChoices = 1:numel(chooseLambda.targetLambdas);
nLambdas = 30; % increase if you're patient
X = trainNEURAL;
ortholambdas = sort([logspace(-1,-7,nLambdas)], 'ascend'); 
% Preallocate
clear loadings cost regularization Hs Ws power
%loadings = [];
%cost  = zeros(length(W1penalties),length(lambdas));
%regularization  = zeros(length(W1penalties),length(lambdas));
%[Hs, Ws] = deal(cell(numel(lambdas),numel(W1penalties)));
%power  = zeros(length(W1penalties),length(lambdas));

[N,T] = size(X);
[OLAM, FITS] = meshgrid(1:numel(ortholambdas), 1:numel(paramChoices));

prog = ProgressBar(2 * numel(OLAM), 'Title', 'Choosing lambda')

dissimilarity_index = 2; 

while dissimilarity_index <= 4
    for use_lam = [1, 0]
        for row = 1:numel(OLAM)
                if size(Ws,1) >= dissimilarity_index && ~isempty(Ws{dissimilarity_index, use_lam+1, row})
                    disp("Continuing")
                    continue
                end
                display(['OrthoH opt: Testing paramset ' num2str(row) '/' num2str(numel(row)*2)])
                % Indices and parameters to describe the computation in this cycle
                i = FITS(row)
                orthoi = OLAM(row);
                W1penality = W1penalties(i);
                if use_lam
                    lambda = targetLambdas(i);
                else
                    lambda = 0;
                end
                % SeqNMF
                [Ws{dissimilarity_index, use_lam+1, row}, Hs{dissimilarity_index, use_lam+1, row}, ...
                 C, loadings(dissimilarity_index, use_lam+1, row,:),...
                 power(dissimilarity_index, use_lam+1, row)] = ...
                    seqNMF_gpu(X, ...
                    'K',single(K),...
                    'L',single(L),...
                    'lambda', lambda,...
                    'lambdaOrthoH', single(ortholambdas(li)),...
                    'maxiter', 150,...
                    'showPlot', 0, ...
                    'lambdaL1H', single(W1penality)); 
                % Compute the cost and reugularization penalty
                [cost(dissimilarity_index, use_lam+1, row),...
                 regularization(dissimilarity_index, use_lam+1, row), ~] =... 
                 helper.get_seqNMF_cost(X,Ws{dissimilarity_index, use_lam+1, row},...
                                        Hs{dissimilarity_index, use_lam+1, row});
                % Save the results
                chooseLambdaOrthoH = struct('loadings', loadings, ...
                                            'power', power, ...
                                            'cost', cost,...
                                            'regularization', regularization, ...
                                            'W1penalties', W1penalties,...
                                            'lambdas', lambdas,...
                                            'orthoi', orthoi,...
                                            'i', 1:numel(OLAM),...
                                            'use_lam', [1,0]);
        end
    end
    mfile.chooseLambdaOrthoH = chooseLambdaOrthoH;
    mfile.chooseLambdaOrthoH_Hs = Hs;
    mfile.chooseLambdaOrthoH_Ws = Ws;
    dissimilarity_index = dissimilarity_index + 1;
end

% Okay now we need to compute dissimilarity
Diss = zeros(size(Hs));
inds = nchoosek(1:dissimilarity_index,2);
for comp = 1:size(inds,1)
    for i = 1:size(Hs,2)
        for j = 1:size(Hs,3)
            disp(comp + "," + i + "," + j);
            Diss(comp,i,j) = ...
            helper.DISSX(Hs{inds(i,1),i,j},...
                         Ws{inds(i,1),i,j},...
                         Hs{inds(i,2),i,j},...
                         Ws{inds(i,2),i,j},'gpu',true);
        end
    end
end

[OLAM, FITS] = ndgrid(1:numel(ortholambdas), 1:numel(paramChoices));

% Visualize dissimilarity surfaces!
muDiss = squeeze(mean(Diss));
% Color use_lam = 0
sz = size(OLAM);
msz = max(sz);
CO1(:,:,1) = zeros(sz); % red
CO1(:,:,2) = ones(sz).*linspace(0.5,0.6,msz); % green
CO1(:,:,3) = ones(sz).*linspace(0,1,msz); % blue
% Color use_lam = 1
CO2(:,:,3) = zeros(sz); % red
CO2(:,:,1) = ones(sz).*linspace(0.5,0.6,msz); % green
CO2(:,:,2) = ones(sz).*linspace(0,1,msz); % blue
fig('LambdaOrthoH-no-smoothing');clf
s1 = surf(ortholambdas(OLAM), paramChoices(FITS), reshape(muDiss(1,:), sz), CO1);
xlabel('\lambda^{\perp}_{H}')
ylabel('\lambda & W1 choice')
%s1 = surf(ortholambdas(OLAM), paramChoices(FITS), reshape(muDiss(1,:), sz), CO1);
set(gca,'xscale','log')
hold on
s2 = surf(ortholambdas(OLAM), paramChoices(FITS), reshape(muDiss(2,:), sz), CO2);
%s2 = surf(ortholambdas(OLAM), paramChoices(FITS), reshape(muDiss(2,:), sz), CO2);
set(gca,'xscale','log')
alpha(s1,0.5);
alpha(s2,0.5);
legend([s1, s2], 'Overall Lambda = false', 'Overall Lambda = true');


% Visualize dissimilarity surfaces!
muDiss = squeeze(mean(Diss));
% Color use_lam = 0
sz = size(OLAM);
msz = max(sz);
CO1(:,:,1) = zeros(sz); % red
CO1(:,:,2) = ones(sz).*linspace(0.5,0.6,msz); % green
CO1(:,:,3) = ones(sz).*linspace(0,1,msz); % blue
% Color use_lam = 1
CO2(:,:,3) = zeros(sz); % red
CO2(:,:,1) = ones(sz).*linspace(0.5,0.6,msz); % green
CO2(:,:,2) = ones(sz).*linspace(0,1,msz); % blue
fig('LambdaOrthoH-smoothing');clf
s1 = surf(ortholambdas(OLAM), paramChoices(FITS), squeeze(smooth3(shiftdim(reshape(muDiss(1,:), sz),-1))), CO1);
xlabel('\lambda^{\perp}_{H}')
ylabel('\lambda & W1 choice')
%s1 = surf(ortholambdas(OLAM), paramChoices(FITS), reshape(muDiss(1,:), sz), CO1);
set(gca,'xscale','log')
hold on
s2 = surf(ortholambdas(OLAM), paramChoices(FITS), squeeze(smooth3(shiftdim(reshape(muDiss(2,:), sz),-1))), CO2);
%s2 = surf(ortholambdas(OLAM), paramChoices(FITS), reshape(muDiss(2,:), sz), CO2);
set(gca,'xscale','log')
alpha(s1,0.5);
alpha(s2,0.5);
legend([s1, s2], 'Overall Lambda = false', 'Overall Lambda = true');

regularizationSurfaces = containers.Map('KeyType', 'int32', 'ValueType', 'any');
costSurfaces = containers.Map('KeyType', 'int32', 'ValueType', 'any');
clear CO1 CO2
[OLAM, FITS] = meshgrid(1:numel(ortholambdas), 1:numel(paramChoices));
smoothit=false;
for use_lam = [0,1]
    % Visualize dissimilarity surfaces!
    sz = size(OLAM);
    msz = max(sz);
    muReg = squeeze(mean(regularization(:,use_lam+1,:)));
    muCost = squeeze(mean(cost(:,use_lam+1,:)));
    muReg = reshape(muReg, sz);
    muCost = reshape(muCost, sz);
    if smoothit
        muReg = squeeze(smooth3(shiftdim(muReg,-1)));
        muCost = squeeze(smooth3(shiftdim(muCost,-1)));
    end
    % Normalize
    muCost=(muCost-min(muCost,[],'all'))./(max(muCost,[],'all')-min(muCost,[],'all'));
    muReg=(muReg-min(muReg,[],'all'))./(max(muReg,[],'all')-min(muReg,[],'all'));
    % Color use_lam = 0
    CO1(:,:,1) = zeros(sz); % red
    CO1(:,:,2) = (ones(sz).*linspace(0.5,0.6,msz)); % green
    CO1(:,:,3) = (ones(sz).*linspace(0,1,msz)); % blue
    % Color use_lam = 1
    CO2(:,:,3) = zeros(sz); % red
    CO2(:,:,1) = (ones(sz).*linspace(0.5,0.6,msz)); % green
    CO2(:,:,2) = (ones(sz).*linspace(0,1,msz)); % blue
    f=fig("CostReg uselam=" + use_lam);figure(f);clf
    cs = reshape(muReg(:), sz);
    costSurfaces(use_lam) = cs;
    s1 = surf(ortholambdas(OLAM), paramChoices(FITS), costSurfaces(use_lam), CO1);
    xlabel('\lambda^{\perp}_{H}')
    ylabel('\lambda & W1 choice')
    %s1 = surf(ortholambdas(OLAM), paramChoices(FITS), reshape(muDiss(1,:), sz), CO1);
    set(gca,'xscale','log')
    hold on
    rs = reshape(muCost(:), sz);
    regularizationSurfaces(use_lam) = rs;
    s2 = surf(ortholambdas(OLAM), paramChoices(FITS), rs, CO2);
    OL = ortholambdas(OLAM);
    PC = paramChoices(FITS);
    %s2 = surf(ortholambdas(OLAM), paramChoices(FITS), reshape(muDiss(2,:), sz), CO2);
    set(gca,'xscale','log')
    alpha(s1,0.5);
    alpha(s2,0.5);
    legend([s1, s2], 'Regularization', 'Cost');
end

% (2) L = R + C 

%% Determine optimal pair by solving for line that parametrizes
% intersection between R and C surfaces, minimizing C!
% (Or alternatively, minimizing C + R

%% First, I have to interpoualte the surface values to get a line through
%resolution_increase = [10, 100];
%regSq = containers.Map('KeyType', 'int32', 'ValueType', 'any');
%costSq = containers.Map('KeyType', 'int32', 'ValueType', 'any');
%method = 'linear'
%for use_lam = [0, 1]
%    OL = ortholambdas(OLAM);
%    PC = paramChoices(FITS);
%    orthoIq      = logspace(ortholambdas(1), ortholambdas(end), numel(ortholambdas)*resolution_increase(1));
%    Iq           = linspace(paramChoices(1), paramChoices(end), numel(paramChoices)*resolution_increase(2));
%    [OIQ, IQ] = meshgrid(orthoIq, Iq);
%    regSq(use_lam) = interp2(OL, PC, double(regularizationSurfaces(use_lam)), OIQ, IQ, method);
%    costSq(use_lam) = interp2(OL, PC, double(costSurfaces(use_lam)), OIQ, IQ, method);
%end
%
%clear CO1 CO2
%for use_lam = [0,1]
%    % Color use_lam = 0
%    sz = size(OIQ);
%    msz = max(sz);
%    CO1(:,:,1) = zeros(sz); % red
%    CO1(:,:,2) = (ones(sz).*linspace(0.5,0.6,msz)'); % green
%    CO1(:,:,3) = (ones(sz).*linspace(0,1,msz)'); % blue
%    % Color use_lam = 1
%    CO2(:,:,3) = zeros(sz); % red
%    CO2(:,:,1) = (ones(sz).*linspace(0.5,0.6,msz)'); % green
%    CO2(:,:,2) = (ones(sz).*linspace(0,1,msz)'); % blue
%    fig("Interp CostReg uselam=" + use_lam);clf
%    cs = costSq(use_lam);
%    s1 = surf(OIQ, IQ, cs, CO1);
%    xlabel('\lambda^{\perp}_{H}')
%    ylabel('\lambda & W1 choice')
%    %s1 = surf(ortholambdas(OLAM), paramChoices(FITS), reshape(muDiss(1,:), sz), CO1);
%    set(gca,'xscale','log')
%    hold on
%    rs = regSq(use_lam);
%    s2 = surf(OIQ, IQ, rs, CO2);
%    OL = ortholambdas(OLAM);
%    PC = paramChoices(FITS);
%    %s2 = surf(ortholambdas(OLAM), paramChoices(FITS), reshape(muDiss(2,:), sz), CO2);
%    set(gca,'xscale','log')
%    alpha(s1,0.5);
%    alpha(s2,0.5);
%    legend([s1, s2], 'Regularization', 'Cost');
%end


% Filtration parameters
windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
% Plot information
I      = chooseLambdaOrthoH.i;
orthoI = chooseLambdaOrthoH.orthoi;
R      = chooseLambdaOrthoH.regularization;
C      = chooseLambdaOrthoH.cost;
% Normalize and smooth
mRs = [prctile(R(:),10), prctile(R(:),90)]; 
mCs = [prctile(C(:),10), prctile(C(:),90)]; 
R = (R-mRs(1))./range(mRs);
C = (C-mCs(1))./range(mCs);
Rs = filtfilt(b,a,double(R)')';
Cs = filtfilt(b,a,double(C)')';
% Plot
figure; clear s
s(1) = surf(I, orthoI, R, 'FaceColor', 'red')
hold on
s(2) = surf(I, orthoI, C, 'FaceColor', 'blue')
legend(s, ["Regularization", "Cost"])
% Plot
figure; clear s
s(1) = surf(I, orthoI, Rs, 'FaceColor', 'red')
hold on
s(2) = surf(I, orthoI, Cs, 'FaceColor', 'blue')
legend(s, ["Regularization", "Cost"])


%                                                                               
% ,---.     |        |   /                  |    |              |        |      
% `---.,---.|---     |__/     ,---.,---.,---|    |    ,---.,-.-.|---.,---|,---.o
%     ||---'|        |  \     ,---||   ||   |    |    ,---|| | ||   ||   |,---| 
% `---'`---'`---'    `   `    `---^`   '`---'    `---'`---^` ' '`---'`---'`---^o
%                                                                               
%                                                                            
%     |     |                   o                             |              
% ,---|,---.|--- ,---.,---.,-.-..,---.,---.    ,---..   .,-.-.|---.,---.,---.
% |   ||---'|    |---'|    | | |||   ||---'    |   ||   || | ||   ||---'|    
% `---'`---'`---'`---'`    ` ' '``   '`---'    `   '`---'` ' '`---'`---'`    
%                                                                            
%                                                                  
%      ,---.         o         ,---.          |                   |
% ,---.|__.     ,---..,---.    |__. ,---.,---.|--- ,---.,---.,---.|
% |   ||        `---.||   |    |    ,---||    |    |   ||    `---. 
% `---'`        `---'``---|    `    `---^`---'`---'`---'`    `---'o
%                     `---'                                        
% %% choose lambda=.005; run multiple times, see number of sig factors
loadings = [];
pvals = []; 
is_significant = []; 
X = trainNEURAL;
nIter = 20; % increase if patient
display('Running seqNMF multiple times for lambda=0.005')

for iteri = 1:nIter
    [W, H, ~,loadings(iteri,:),power]= seqNMF(X,'K',K,'L',Lneural,...
            'lambdaL1H', .1, 'lambda', targetLambda, 'maxiter', 100, 'showPlot', 0); 
    p = .05;
    [pvals(iteri,:),is_significant(iteri,:)] = test_significance(testNEURAL,W);
    W = W(:,is_significant(iteri,:)==1,:); 
    H = H(is_significant(iteri,:)==1,:); 
    [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
    indSort = hybrid(:,3);
    tstart = 300; 
    clf; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 0,...
                 trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end))
    display(['seqNMF run ' num2str(iteri) '/' num2str(nIter)])
end
figure; hold on
h = histogram(sum(is_significant,2), 'edgecolor', 'w', 'facecolor', .7*[1 1 1]); 
h.BinCounts = h.BinCounts/sum(h.BinCounts)*100; 
xlim([0 10]); 
xlabel('# significant factors')
ylabel('% seqNMF runs')

%% Plot factor-triggered song examples and rastors
addpath(genpath('misc_elm')); 
figure; HTriggeredSpec(H,trainSONG,VIDEOfs,SONGfs,Lsong); 
%figure; HTriggeredRaster(H,trainNEURAL(indSort,:),Lneural);

%  ____           _                            _  __
% |  _ \ ___  ___| |__   ___   ___  ___  ___  | |/ /
% | |_) / _ \/ __| '_ \ / _ \ / _ \/ __|/ _ \ | ' / 
% |  _ <  __/ (__| | | | (_) | (_) \__ \  __/ | . \ 
% |_| \_\___|\___|_| |_|\___/ \___/|___/\___| |_|\_\
%                                                   
%% Procedure for choosing K
% Given that we optomized lambda, let's choose K 1 more time!
tic
Ws = {};
Hs = {};
numfits = 3; %number of fits to compare
for k = 1:12
    display(sprintf('running seqNMF with K = %i',k))
    for ii = 1:numfits
        [Ws{ii,k},Hs{ii,k}] = seqNMF_gpu(X,'K',k, 'L', L,'lambda', 0,'maxiter',100,'showplot',0); 
        % note that max iter set low (30iter) for speed in demo (not recommended in practice)
    end
    inds = nchoosek(1:numfits,2);
    for i = 1:size(inds,1) % consider using parfor for larger numfits
            Diss(i,k) = helper.DISSX(Hs{inds(i,1),k},Ws{inds(i,1),k},Hs{inds(i,2),k},Ws{inds(i,2),k});
    end
    
end
toc
%% Plot Diss and choose K with the minimum average diss.
figure,
plot(1:10,Diss,'ko'), hold on
h1 = plot(1:10,median(Diss,1),'k-','linewidth',2);
h2 = plot([3,3],[0,0.5],'r--');
legend([h1 h2], {'median Diss','true K'})
xlabel('K')
ylabel('Diss')

display('Testing significance of factors on held-out data')
[pvals,is_significant] = test_significance(testNEURAL,W);

W = W(:,is_significant,:); 
H = H(is_significant,:); 

%% And maybe just once more examine significant factors with lambdaOrtho, even though it
% reduces simiarity of the Hs and Ws found.

%% Graveyard
%
% Dissimilarity LambdaOrthoH section
% ----------------------------------
%Code to fix improper save structure
%chooseLambdaOrthoH = mfile.chooseLambdaOrthoH;
%Ws = reshape({chooseLambdaOrthoH.Ws}, size(chooseLambdaOrthoH));
%Hs = reshape({chooseLambdaOrthoH.Hs}, size(chooseLambdaOrthoH));
%chooseLambdaOrthoH = chooseLambdaOrthoH(end);
%chooseLambdaOrthoH.i           = 1:size(chooseLambdaOrthoH.regularization,2);
%chooseLambdaOrthoH.orthoi      = 1:size(chooseLambdaOrthoH.regularization,1);
%chooseLambdaOrthoH.W1penalties = W1penalties(FITS(:));
%chooseLambdaOrthoH.lambdas     = targetLambdas(FITS(:));
%mfile.chooseLambdaOrthoH = chooseLambdaOrthoH;
%
%%Code to fix improper save structure
%% (Prepping for adding dissimilarity dimension)
%chooseLambdaOrthoH = mfile.chooseLambdaOrthoH;
%Ws = reshape({chooseLambdaOrthoH.Ws}, size(chooseLambdaOrthoH));
%Hs = reshape({chooseLambdaOrthoH.Hs}, size(chooseLambdaOrthoH));
%Ws = shiftdim(Ws{1},-1);
%Hs = shiftdim(Hs{1},-1);
%chooseLambdaOrthoH = chooseLambdaOrthoH(end);
%chooseLambdaOrthoH.i           = 1:size(chooseLambdaOrthoH.regularization,2);
%chooseLambdaOrthoH.orthoi      = 1:size(chooseLambdaOrthoH.regularization,1);
%chooseLambdaOrthoH.dissim      = 1;
%chooseLambdaOrthoH.W1penalties = W1penalties(FITS(:));
%chooseLambdaOrthoH.lambdas     = targetLambdas(FITS(:));
%chooseLambdaOrthoH.regularization = shiftdim(chooseLambdaOrthoH.regularization,-1);
%chooseLambdaOrthoH.cost           = shiftdim(chooseLambdaOrthoH.cost,-1);
%chooseLambdaOrthoH.loadings       = shiftdim(chooseLambdaOrthoH.loadings,-1);
%chooseLambdaOrthoH.power          = shiftdim(chooseLambdaOrthoH.power,-1);
%for field = fieldnames(chooseLambdaOrthoH)'
%    field{1}
%    assign(field{1}, chooseLambdaOrthoH.(field{1}));
%end
%dmfile.chooseLambdaOrthoH = chooseLambdaOrthoH;

