% Things to do
% -------------
% 1. optomize the lambda
% 2. optomize L? (lower priority, seems to remove patterns on its own -- possibly speeds)
% 3. Another animal
%   a. Correlate
%   b. W(animal_1) = W(animal_2) in terms of loading
% 4. Spectral frequencies of a W?

% Boilerplate
% -----------
seqnmf_initialize; % Boilerplate across analyses
% Lower number of iterations for convergence
for i = 1:numel(seqnmf_kws)
    f = find(cellfun(@(x) isequal(x, 'maxiter'), seqnmf_kws{i}));
    seqnmf_kws{i}{f+1} = 40; % According to cost function progression, more than enough iterations to asymptotically minimize
end
folder = [folder filesep 'crossvalidate'];
disp(['Folder = ' folder]);
nReps = 100;

assert(~(ismember('animal_overall', nmfs) && ismember('overall', nmfs)))

 % Add animal data to path
% ------------------------
for animal = animal_list
    animalToPath(animal{1});
end

disp("Starting " + string(seqnmf_arg_names(argNum)))
% Iterate each animal and run seqnmf logic
% ----------------------------------------
P = []; % for aggreagting pos
C = []; % for aggreagting pos
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
        L =  round(10/data.params.movingwin(2));

        if ismember('epoch', nmfs)
            data  = ry_run_seqnmf(data, fields, 'seqnmf_kws', {seqnmf_kws{argNum}{:}, 'L', L, 'checkpoint', 200}, 'fieldpack_kws', fieldpack_kws);
        end
    end
end

file = @(x) sprintf('~/Data/deltacoherence/%s/%s/master_seqnmf_rep=%d.mat', folder, seqnmf_arg_names(argNum), x);
file = @(x) sprintf('~/Data/deltacoherence/%s/%s/master_seqnmf_rep=%d.mat', folder, seqnmf_arg_names(argNum), x);
% Aggregate
% ---------
data.t = cat(2, C.t);
data.f = cat(1, C.f);
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

RmseTrain = zeros(1, nReps);
RmseTest = zeros(1, nReps);

% Examine
% -------
sfigure; 
clf
%kPlot = 1:nReps;
b1= bar(1, mean(RmseTrain), 'FaceColor', 'b')
jitter = - 0.5 + 1 * rand(size(RmseTrain));
scatter(1*ones(size(RmseTrain)) + jitter, RmseTrain(:), 'k', 'markerfacecolor', 'flat'); 
hold on; 
b2= bar(2, mean(RmseTrain), 'FaceColor', 'b')
scatter(2*ones(size(RmseTest)) + jitter, RmseTest(:), 'b', 'markerfacecolor', 'flat'); 
plot(mean(RmseTrain,2), 'r')
plot(mean(RmseTest,2), 'b')
xlabel('K'); ylabel('RMSE')
legend('Train', 'Test', 'location', 'northwest')
drawnow; shg
%linkdata on;


% Cross-validation analysis here
% ------------------------------
sfigure;
data.data = seqnmf_packfields(data, fields, fieldpack_kws{:});
X = data.data';
tic; tic;
for repi = 1:nReps
    display(['Cross validation on masked test set; rep ' num2str(repi)])
    rng('shuffle')
    M = rand(size(X,1),size(X,2)) > .05; % create masking matrix (0's are test set, not used for fit)
    [W, H] = seqnmf_factorize(X,'M', M, seqnmf_kws{argNum}{:});
    Xhat = helper.reconstruct(W,H); 
    RmseTrain(1,repi) =   sqrt(sum(M(:).*(X(:)-Xhat(:)).^2)  ./ sum(M(:)));
    RmseTest( 1, repi) =  sqrt(sum((~M(:)).*(X(:)-Xhat(:)).^2)./sum(~M(:)));
    PowerTrain(1,repi) =  (sum((M(:).*X(:)).^2)-sum((M(:).*(X(:)-Xhat(:))).^2))/sum((M(:).*X(:)).^2);  % fraction power explained by whole reconstruction
    PowerTest( 1, repi) = (sum((~M(:).*X(:)).^2)-sum((~M(:).*(X(:)-Xhat(:))).^2))/sum((~M(:).*X(:)).^2);  % fraction power explained by whole reconstruction
    if mod(repi, 10) == 0
        data.W = W;
        data.H = H;
        save(file(repi), '-struct', 'data')
    end
    toc
end
toc

% Plot final result
% -----------------
colors = cmocean('curl',6);
ctest = colors(2,:);
ctrain = colors(4,:);
sfigure;
bincount = 6;
[~, train_edges] = histcounts(RmseTrain, bincount);
[~, test_edges]  = histcounts(RmseTest, bincount);
% ---------
subplot(2,2,1);
histogram(RmseTrain, train_edges, 'FaceColor', ctrain, 'FaceAlpha', 0.25);
hold on;
v= vline(median(RmseTrain)); v.Color = ctrain;
histogram(RmseTest, test_edges, 'FaceColor', ctest, 'FaceAlpha', 0.25)
v= vline(median(RmseTest)); v.Color = ctest;
l_extent = [median(RmseTrain), median(RmseTest)];
l = line(l_extent, [25 25], 'Color','white', 'LineStyle', '--')
difference_text = sprintf('diff = %1.2e', diff(l_extent));
text(l_extent(1) + 0.25*range(l_extent), 26, difference_text, 'fontsize', 20)
xlabel('Mean squared error over cross-validations')
legend('Train','Test')
% ---------
RTr = 100 * RmseTrain / 3;
RTe = 100 * RmseTest / 3;
test_edges = 100 * test_edges / 3;
train_edges = 100 * train_edges / 3;
subplot(2,2,3);
histogram(RTr, train_edges, 'FaceColor', ctrain, 'FaceAlpha', 0.25);
hold on;
v= vline(median(RTr)); v.Color = ctrain;
histogram(RTe, test_edges, 'FaceColor', ctest, 'FaceAlpha', 0.25)
v= vline(median(RTe)); v.Color = ctest;
l_extent = [median(RTr), median(RTe)];
l = line(l_extent, [25 25], 'Color','white', 'LineStyle', '--')
difference_text = sprintf('diff = %1.2e', diff(l_extent));
text(l_extent(1) + 0.25*range(l_extent), 26, difference_text, 'fontsize', 20)
xlabel('Percent RMSE of total value range (cross-val)')
legend('Train','Test')
% ---------
subplot(2,2,[2 4]); cla
jitter = -0.5 + rand(size(PowerTrain));
s = scatter(ones(size(PowerTrain))*1 + jitter, PowerTrain)
set(s, 'MarkerFaceColor', ctrain, 'MarkerFaceAlpha', 0.50, 'MarkerEdgeColor', ctrain, 'MarkerEdgeAlpha', 0.5);
s= scatter(ones(size(PowerTrain))*1 + jitter, PowerTest)
set(s, 'MarkerFaceColor', ctest, 'MarkerFaceAlpha', 0.50, 'MarkerEdgeColor', ctrain, 'MarkerEdgeAlpha', 0.5);
legend('Power_{train}', 'Power_{test}')
ylabel('fraction \sigma^2 explained (power)');

save(sprintf('~/Projects/deltacoherence/results/%s/%s/master_crossval.mat', folder, seqnmf_arg_names(argNum)), 'RmseTest','RmseTest', 'PowerTest', 'PowerTrain');

saveas(gcf, sprintf('~/Projects/deltacoherence/results/%s/%s/master_crossval.svg', folder, seqnmf_arg_names(argNum)));
saveas(gcf, sprintf('~/Projects/deltacoherence/results/%s/%s/master_crossval.fig', folder, seqnmf_arg_names(argNum)));
figure; imagesc(M');set(gca,'ydir','normal'); colormap('gray'); title("Example M"); saveas(gcf,sprintf('~/Projects/deltacoherence/results/%s/%s/example_m.fig', folder, seqnmf_arg_names(argNum)));
