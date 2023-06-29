
seqnmf_initialize

% Add animal data to path
% ------------------------
for animal = animal_list
    animalToPath(animal{1});
end

% Iterate each animal and run seqnmf logic
% ----------------------------------------
searchFolder = '~/Data/deltacoherence/'
pushd(searchFolder)
files = subdir('*seqnmf.mat')
figure_count = 1;
for file = files'

    folder =  file.folder;
    tic
    data = load(file.name);
    disp(['Load required ' num2str(toc) ' seoncds']);

    [path, ~] = fileparts(file.name);
    [~, foldername] = fileparts(path);

    % If data does not have fields
    if ~isfield(data,'fields')
        warning('File did not contain a fields set. Creating one.')
        fields = {};
        maxnum = max(structfun(@numel, data));
        for field = fieldnames(data)'
            if numel(data.(field{1}) == maxnum
                fields{end+1} = field;
            end
        end
        data.fields = fields;
    end

    %% REPEATED PROPERTIES
    % foldername : this is the key string sequence distinguishing these analyses from eachother
    % figure_count : unique number for that figure

    %%% Data
    %data_plot = figure(100 + figure_count); clf
    %fstring = string(F);
    %file = sprintf('~/Data/deltacoherence/%s/%s/data', folder, seqnmf_arg_names(argNum));
    %if isequal(epoch_type,'run')
    %    yscale = size(master.data,1);
    %    imagesc(master.time, 1:yscale, master.data'); set(gca,'ydir','normal')
    %    cmocean('thermal')
    %    hold on
    %    plot(master.time, yscale*(master.trajdist + 1)/2, 'k', 'LineStyle', ':', 'LineWidth', 3)
    %    plot(master.time, smoothdata(yscale*master.vel./prctile(master.vel, 90)), 'LineStyle','--', 'Color', 'white')
    %else
    %    yscale = size(data.data,1);
    %    imagesc(master.time, 1:yscale, master.data'); set(gca,'ydir','normal')
    %    cmocean('thermal')
    %    hold on
    %end
    %xlim([6.3519    6.5982]*1000);
    %suptitle(sprintf('%s -- Timescale: %d, fields: %s', foldername, timescale, fstring.join('-')));
    %saveas(data_plot, sprintf('%s.%s', file, 'svg'))
    %try
    %    saveas(data_plot, sprintf('%s.%s', file, 'fig'))
    %catch ME
    %    save(sprintf('%s.%s', file, 'fig.mat'), 'data_plot', '-v7.3')
    %end

    %%% Learned patterns
    %pattern_plot = figure(1001 + figure_count);
    %fstring = string(F);
    %file = sprintf('~/Data/deltacoherence/%s/%s/master_seqnmf', folder, seqnmf_arg_names(argNum));
    %%ry_WHPlot(W, H, data, 'plotAll', 0, 'ExtraMatrixToPlot', WH, 'plotWflat', 0);
    %ry_WHPlot(data.W, data.H, data.data, 'plotAll', 0, 'ExtraMatrixToPlot', data.WH, 'plotWflat', 0);
    %initial = sprintf('%s, Timescale: %d, fields: %s', foldername, timescale, fstring.join('-'));
    %suptitle(sprintf('%s, cost=%2.4f, power=%2.4f\nlambda=0.0009', initial, animal, mean(cost(end)), power))
    %
    %saveas(pattern_plot, sprintf('%s.%s', file, 'svg'))
    %saveas(pattern_plot, sprintf('%s.%s', file, 'fig'))
    
    %%% Ws
    % W_plot
    seqnmf_plotdata(data, 'fields', data.fields), 'plotType', 'wplot')

    figure_count = figure_count + 1
end

