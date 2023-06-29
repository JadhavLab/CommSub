rootFolder = '/Volumes/Data/deltacoherence/';
skipMaster = true;

% The name of the root folder
rootedfolder = @(x) fullfile(rootFolder, x);

% Replots all of the results from SeqNMF
InitializeParamsets;
% Iterates over each of the parameter sets and replots
% the data
diary(['~/Projects/deltacoherence/results/diary_replot_' date])
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
    if isempty(folderfull); 
        disp([fullfile(rootFolder,folder) ' is empty'])
        continue; 
    else
        data = M.data;
        data.WH = helper.reconstruct(data.W, data.H)
    end

    % It's hammer time! - MC Plotter
    % ----------------------------
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

    % LEARNED PATTERN
    % ---------------
    pattern_plot = fig(fullfile(folderfull, 'learned pattern'));
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    fstring = string(data.fields);
    file = fullfile(folderfull, seqStyle, 'master_seqnmf');
    %ry_WHPlot(W, H, data, 'plotAll', 0, 'ExtraMatrixToPlot', WH, 'plotWflat', 0);
    [axP, axR, axisCenters] = labeler.wholeLabel(data,data.fields, 'downsample_specOnly', 20);
    [~,indplot] = ry_WHPlot(data.W, data.H, data.data, 'wwidth', 0.35, 'plotAll', 0, 'ExtraMatrixToPlot', data.WH, 'plotWflat', 0, 'Fs', 1/median(diff(data.t)),...
        'yaxis_centers', axisCenters, 'yaxis_axP', axP);
    suptitle(sprintf('ArgNum: %s, Timescale: %d, fields: %s', seqStyle, timescale, fstring.join('-')));
    saveas(pattern_plot, sprintf('%s.%s', file, 'svg'))
    saveas(pattern_plot, sprintf('%s.%s', file, 'fig'))

    % Loadings
    loading_plot = fig(char(fullfile('~/Projects/deltacoherence/results', folder, seqStyle, 'loadings')));
    plot(1:K, data.loadings, '*-');

    % LEARNED PATTERN raw only wo wpli
    % ---------------
    pattern_plot = fig(fullfile(folderfull, 'learned pattern wo raw wo wpli'));
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %fstring = string(F);
    file = fullfile(folderfull, seqStyle, 'learnedpatternWOraw');
    %ry_WHPlot(W, H, data, 'plotAll', 0, 'ExtraMatrixToPlot', WH, 'plotWflat', 0);
    %ry_WHPlot(data.W, data.H, data.data, 'wwidth', 0.35, 'plotAll', 0, 'plotWflat', 0, 'Fs', 1/median(diff(data.t)),...
    %    'yaxis_centers', axisCenters, 'yaxis_axP', axP, 'imgaussfiltscale', 0);
    %meanWgroup = {sieve.field(data,'wpli')};
    meanWgroup = {};
    remove = ismember(data.fields, {'wpli'});
    remove = data.fields(remove);
    removeWgroup = arrayfun(@(x) sieve.field(data,'wpli'), remove, 'UniformOutput', false);
    axP=replace(axP, '$','');
    [axP, axR, axisCenters, fieldPts] = labeler.wholeLabel(data, data.fields(~ismember(data.fields,'wpli')));
    colorWgroups = arrayfun(@(x) axisCenters(fieldPts == x), ["velocity", "trajdist"], 'UniformOutput', false);
    [axP, axR, axisCenters, fieldPts] = labeler.wholeLabel(data, data.fields(~ismember(data.fields,'wpli')), 'downsample_specOnly', 20);
    axP=replace(axP, '$','');
    colorWgroups = arrayfun(@(x) axisCenters(fieldPts == x), ["velocity", "trajdist"], 'UniformOutput', false);
    colorWmaps = {'haline','phase'};
    [imCell, indplot] = ry_WHPlot(data.W, data.H, data.data, 'wwidth', 0.35, 'plotAll', 0, 'plotWflat', 1, 'Fs', 1/median(diff(data.t)),...
        'yaxis_centers', axisCenters, 'yaxis_axP', axP, 'imgaussfiltscale', 0, 'meanWgroup', meanWgroup, 'removeWgroup', removeWgroup, 'colorWgroups', colorWgroups, 'colorWmaps', colorWmaps);
        fstring = string(data.fields);
    suptitle(sprintf('ArgNum: %s, Timescale: %d, fields: %s', seqStyle, data.timescale, fstring.join('-')));
    saveas(pattern_plot, sprintf('%s_nowpli.%s', file, 'svg'))
    saveas(pattern_plot, sprintf('%s_nowpli.%s', file, 'fig'))

    % LEARNED PATTERN raw only wo wpli pfc
    % ---------------
    pattern_plot = fig(fullfile(folderfull, 'learned pattern wo raw wo wpli pfc'));
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %fstring = string(F);
    file = fullfile(folderfull, seqStyle, 'learnedpatternWOraw');
    %ry_WHPlot(W, H, data, 'plotAll', 0, 'ExtraMatrixToPlot', WH, 'plotWflat', 0);
    %ry_WHPlot(data.W, data.H, data.data, 'wwidth', 0.35, 'plotAll', 0, 'plotWflat', 0, 'Fs', 1/median(diff(data.t)),...
    %    'yaxis_centers', axisCenters, 'yaxis_axP', axP, 'imgaussfiltscale', 0);
    meanWgroup = {};
    remove = ismember(data.fields, {'wpli','S2'});
    remove = data.fields(remove);
    removeWgroup = arrayfun(@(x) sieve.field(data,x), remove, 'UniformOutput', false);
    [axP, axR, axisCenters, fieldPts] = labeler.wholeLabel(data, data.fields(~ismember(data.fields,{'wpli','S2'})));
    axP=replace(axP, '$','');
    colorWgroups = arrayfun(@(x) axisCenters(fieldPts == x), ["velocity", "trajdist"], 'UniformOutput', false);
    [axP, axR, axisCenters, fieldPts] = labeler.wholeLabel(data, data.fields(~ismember(data.fields,{'wpli','S2'})), 'downsample_specOnly', 20);
    axP=replace(axP, '$','');
    colorWgroups = arrayfun(@(x) axisCenters(fieldPts == x), ["velocity", "trajdist"], 'UniformOutput', false);
    colorWmaps = {'haline','phase'};
    ry_WHPlot(data.W, data.H, data.data, 'wwidth', 0.35, 'plotAll', 0, 'plotWflat', 1, 'Fs', 1/median(diff(data.t)),...
        'yaxis_centers', axisCenters, 'yaxis_axP', axP, 'imgaussfiltscale', 0, 'meanWgroup', meanWgroup, 'removeWgroup', removeWgroup, 'colorWgroups', colorWgroups, 'colorWmaps', colorWmaps, 'indplot', indplot);
        fstring = string(data.fields);
    suptitle({sprintf('ArgNum: %s, Timescale: %d, fields: %s', seqStyle, data.timescale, fstring.join('-')),...
    sprintf('Power: %2.1f%%, Loadings: %s',data.power*100, sprintf('%0.2f ', data.loadings))});
    saveas(pattern_plot, sprintf('%s_nopfc.%s', file, 'svg'))
    saveas(pattern_plot, sprintf('%s_nopfc.%s', file, 'fig'))

    % LEARNED PATTERN raw only wo wpli hpc
    % ---------------
    pattern_plot = fig(fullfile(folderfull, 'learned pattern wo raw wo wpli hpc'));
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %fstring = string(F);
    file = fullfile(folderfull, seqStyle, 'learnedpatternWOraw');
    %ry_WHPlot(W, H, data, 'plotAll', 0, 'ExtraMatrixToPlot', WH, 'plotWflat', 0);
    %ry_WHPlot(data.W, data.H, data.data, 'wwidth', 0.35, 'plotAll', 0, 'plotWflat', 0, 'Fs', 1/median(diff(data.t)),...
    %    'yaxis_centers', axisCenters, 'yaxis_axP', axP, 'imgaussfiltscale', 0);
    meanWgroup = {};
    remove = ismember(data.fields, {'wpli','S1'});
    remove = data.fields(remove);
    removeWgroup = arrayfun(@(x) sieve.field(data,x), remove, 'UniformOutput', false);
    [axP, axR, axisCenters, fieldPts] = labeler.wholeLabel(data, data.fields(~ismember(data.fields,{'wpli','S1'})));
    axP=replace(axP, '$','');
    colorWgroups = arrayfun(@(x) axisCenters(fieldPts == x), ["velocity", "trajdist"], 'UniformOutput', false);
    [axP, axR, axisCenters, fieldPts] = labeler.wholeLabel(data, data.fields(~ismember(data.fields,{'wpli','S1'})), 'downsample_specOnly', 20);
    axP=replace(axP, '$','');
    colorWgroups = arrayfun(@(x) axisCenters(fieldPts == x), ["velocity", "trajdist"], 'UniformOutput', false);
    colorWmaps = {'haline','phase'};
    ry_WHPlot(data.W, data.H, data.data, 'wwidth', 0.35, 'plotAll', 0, 'plotWflat', 1, 'Fs', 1/median(diff(data.t)),...
        'yaxis_centers', axisCenters, 'yaxis_axP', axP, 'imgaussfiltscale', 0, 'meanWgroup', meanWgroup, 'removeWgroup', removeWgroup, 'colorWgroups', colorWgroups, 'colorWmaps', colorWmaps, 'indplot', indplot);
        fstring = string(data.fields);
    suptitle({sprintf('ArgNum: %s, Timescale: %d, fields: %s', seqStyle, data.timescale, fstring.join('-')),...
    sprintf('Power: %2.1f%%, Loadings: %s',data.power*100, sprintf('%0.2f ', data.loadings))});
    saveas(pattern_plot, sprintf('%s_nohpc.%s', file, 'svg'))
    saveas(pattern_plot, sprintf('%s_nohpc.%s', file, 'fig'))

    % LEARNED PATTERN raw only wo wpli hpc pfc
    % ---------------
    pattern_plot = fig(fullfile(folderfull, 'learned pattern wo raw wo wpli hpc pfc'));clf
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %fstring = string(F);
    file = fullfile(folderfull, seqStyle, 'learnedpatternWOraw');
    %ry_WHPlot(W, H, data, 'plotAll', 0, 'ExtraMatrixToPlot', WH, 'plotWflat', 0);
    %ry_WHPlot(data.W, data.H, data.data, 'wwidth', 0.35, 'plotAll', 0, 'plotWflat', 0, 'Fs', 1/median(diff(data.t)),...
    %    'yaxis_centers', axisCenters, 'yaxis_axP', axP, 'imgaussfiltscale', 0);
    meanWgroup = {};
    remove = ismember(data.fields, {'wpli','S1','S2'});
    remove = data.fields(remove);
    removeWgroup = arrayfun(@(x) sieve.field(data,x), remove, 'UniformOutput', false);
    [axP, axR, axisCenters, fieldPts] = labeler.wholeLabel(data, data.fields(~ismember(data.fields,{'wpli','S1','S2'})));
    axP=replace(axP, '$','');
    colorWgroups = arrayfun(@(x) axisCenters(fieldPts == x), ["velocity", "trajdist"], 'UniformOutput', false);
    [axP, axR, axisCenters, fieldPts] = labeler.wholeLabel(data, data.fields(~ismember(data.fields,{'wpli','S1','S2'})), 'downsample_specOnly', 20);
    axP=replace(axP, '$','');
    colorWgroups = arrayfun(@(x) axisCenters(fieldPts == x), ["velocity", "trajdist"], 'UniformOutput', false);
    colorWmaps = {'haline','phase'};
    ry_WHPlot(data.W, data.H, data.data, 'wwidth', 0.35, 'plotAll', 0, 'plotWflat', 1, 'Fs', 1/median(diff(data.t)),...
        'yaxis_centers', axisCenters, 'yaxis_axP', axP, 'imgaussfiltscale', 0, 'meanWgroup', meanWgroup, 'removeWgroup', removeWgroup, 'colorWgroups', colorWgroups, 'colorWmaps', colorWmaps, 'indplot', indplot);
    fstring = string(data.fields);
    suptitle({sprintf('ArgNum: %s, Timescale: %d, fields: %s', seqStyle, data.timescale, fstring.join('-')),...
    sprintf('Power: %2.1f%%, Loadings: %s',data.power*100, sprintf('%0.2f ', data.loadings))});
    saveas(pattern_plot, sprintf('%s_nohpcpfc.%s', file, 'svg'))
    saveas(pattern_plot, sprintf('%s_nohpcpfc.%s', file, 'fig'))
end
