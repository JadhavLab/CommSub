% ,--.      |    |         
% |   |,---.|    |--- ,---.
% |   ||---'|    |    ,---|
% `--' `---'`---'`---'`---^
%           |                                  ,---.
% ,---.,---.|---.,---.,---.,---.,---.,---.,---.  ,-'
% |    |   ||   ||---'|    |---'|   ||    |---'  |  
% `---'`---'`   '`---'`    `---'`   '`---'`---'  o  
%
% Auth: Ryan Y
%
% If it's real, it should show up in several animals in
% the final epochs.



% Parameters
% ----------
bands               = [6 12; 0.5 4]; % 1- theta, 2 - delta
lowhigh_defintions = [0.2, 0.8]; % percentile definitions for high and low respectively
animal_list         = {'JS12','JS13','JS14','JS15','JS21','ER1','KL8'};
fields              = {'S1','S2','wpli'};
CA1 = 1; PFC = 2; CA1PFC = 3;
theta = 1; delta = 2;
high = 2; low = 1;
frequency = 2; time = 1;
samprate = 0.1;
win = [10, 10]; % 10 seconds before and after
plot_range = [0, 40];



% Add animal data to path
% ------------------------
for animal = animal_list
    animalToPath(animal{1});
end

% Let's carray out this process for every animal
for animal = animal_list

    %% SECTION 1: PERFORM COMPUTATIONS

    % Determine the days and epochs for the animal
    animaldef_ = animaldef(animal{1});
    task = loaddatastruct(animaldef_{2:3}, 'task');
    search_results = cellfetch(task, 'environment');
    dayepochs = search_results.index;

    % Load up the cgramc
    for dayepoch = dayepochs'

        % Load
        % ----
        % load spectral info
        cgramc_file = [animal{1} 'cgramc' sprintf('-%02d-%02d.mat', dayepoch)];
        if ~exist(cgramc_file, 'file')
            continue
        end
        load(cgramc_file)
        data = ffend(cgramc);
        % load position info
        load([animal{1} 'pos' sprintf('%02d', dayepoch(1))])
        day = dayepoch(1);
        epoch = dayepoch(2);
        P = pos{day}{epoch};

        % Throw out noise periods
        % -----------------------
        threshold = 1.5; % standard deviations
        noise = all(zscore( data.S1, 1, 1) > threshold, frequency) ...
        & all(zscore(data.S2, 1, 1) > threshold, frequency);
        data.t(noise) = 0;
        for field = fields
            data.(field{1})(noise,:) = 0;
        end

        %% %%%%%%% Determine certain components %%%%%%%%%%%%%

        % Upfront necessaries : 
        % ----------------------
        % - raw{ spec_hpc | spec_pfc | wpli }{ theta | delta } : 
        %           records the raw values for the datatypes in theta and delta
        % - avg{ spec_hpc | spec_pfc | wpli }{ theta | delta } : 
        %           records bootstrapped averages
        % - event{ theta | delta }{ high | low } : 
        %           records the high and low of these bands in pfc
        
        % (A) Raw
        % raw -- field -- band
        % -------
        field_cnt = 0;
        for field = fields
            band_cnt = 0;
            field_cnt = field_cnt + 1;
            for band = bands'
                band_cnt = band_cnt + 1;
                datafield = data.(field{1});
                filter = data.f >= band(1) ...
                       & data.f < band(2);
                raw{field_cnt}{band_cnt} = mean(datafield(:, filter), frequency);
            end
        end

        % (B) Event
        % event (just pfc spectra right now) --- band -- high/low
        % ----------
        band_cnt = 0;
        for band = bands'
            band_cnt = band_cnt + 1;

            for condition = [low, high]

                if condition == low % low
                    Q = quantile( raw{PFC}{band_cnt}, lowhigh_defintions(condition));
                    quantile_value{band_cnt}{condition} = Q;
                    event{band_cnt}{condition} = ...
                        raw{2}{band_cnt} < Q;
                elseif condition == high % high
                    Q = quantile( raw{PFC}{band_cnt}, lowhigh_defintions(condition));
                    quantile_value{band_cnt}{condition} = Q;
                    event{band_cnt}{condition} = ...
                        raw{2}{band_cnt} > Q;
                end
            end
        end

        % (C) Determine bootstrapped averages
        % averages -- field --- band -- high/low
        % -----------------------------------
        for field_cnt = 1:numel(fields)
            for band = 1:size(bands,1)
                datafield = raw{field_cnt}{band};
                for condition_cnt = [low, high]

                    opt.UseParallel=true;
                    times = event{delta}{condition_cnt};
                    average{field_cnt}{band}{condition_cnt} = ...
                        bootstrp(1000, @mean, datafield(times)',...
                            'Options', opt);
                end
            end
        end

        % (C) Get velocity data timed to spectral
        % ---------------------------------------
        indices  = lookup(data.t, P.data(:,1));
        velocity = P.data(indices, 5);

        % (D) Triggered windows 
        % ---------------------
        for band = 1:size(bands, 1)
        
            [peaks, peak_locations ] = findpeaks(raw{PFC}{band})
            for condition_cnt = 1:numel(lowhigh_defintions)
                if condition_cnt == low
                    locs_of_interest = peak_locations(peaks < quantile_value{band}{low});
                else
                    locs_of_interest = peak_locations(peaks > quantile_value{band}{high});
                end

                % Obtain windows
                for field_cnt = 1:numel(fields)
                    datafield = raw{field_cnt}{band};

                    % Compute indices of windows
                    startstops = ...
                    [locs_of_interest - round(win(1)/samprate), ...
                    locs_of_interest + round(win(2)/samprate)];

                    % Determine bad windows (those extending before or after data)
                    % and delete them
                    windows_before_data_collection  = startstops(:,1) < 1 ;
                    windows_after_data_collection = startstops(:,2) > numel(datafield);
                    startstops( windows_after_data_collection | windows_before_data_collection, :) = [];
                    
                    nWindows = size(startstops,1);
                    winsize = round(win(1)/samprate) + round(win(2)/samprate) + 1;

                    % Collect windows of data around events
                    transform = @(x) zscore(x,1,1);
                    datafield = transform(datafield);
                    windows{field_cnt}{band} = zeros(nWindows, winsize);
                    for window = 1:nWindows
                        windows{field_cnt}{band}(window, :) = ...
                            datafield(startstops(window,1):startstops(window,2));
                    end


                end
            end
        end

        stat_tab = table('RowNames', {'animal','day','epoch','stat','val','area','band','cond'});
        add_row = @(stat,val,area,band,cond) [stat_tab;{animal,day,epoch,stat,val,area,band,cond}];
        % add_mean ;
        % add_mean_cl ;

        %% SECTION 2: PLOT TIME!
        % ---------------------------------------------
        % (1) Determine periods of high and low delta
        figure;
        clf;
        ax = subplot(2,2,1:2);
        plotevents = true;
        title(['Session' sprintf('-%s-%d-%d',animal{1}, dayepoch)])
        F = data.f > plot_range(1) & data.f < plot_range(2);
        data.f = data.f(F);
        norm_f  = @(x) ((x / range(quantile(x(:),[0, 0.8])))-0.25)*range(data.f)/2 + min(data.f);
        for field = fields
            data.(field{1}) = data.(field{1})(:,F);
        end
        nestable on;
        ax_sub = gobjects(3,1);
        for i = 1:3
            spectral_data = real(data.(fields{i}));
            ax_sub(i) =  nestplot(3,1,i);
            fields{i}
            if fields{i}(1) == 'w'
                imagesc(ax_sub(i), data.t, data.f, spectral_data', [-0.2, 0.5]);
            elseif fields{i}(1) == 'S'
                %spectral_data = log10(spectral_data+10);
                spectral_data = zscore(spectral_data,1,1);
                imagesc(ax_sub(i), data.t, data.f, spectral_data', [-1, 2]);
                %set(gca,'clim',[-1,2])
            end
            colorbar
            set(gca,'ydir','normal')
            hold on;
            if isequal(fields{i}, 'S1')
                ylabel(ax_sub(i),'CA1')
            elseif isequal(fields{i}, 'S2')
                ylabel(ax_sub(i),'PFC')
            elseif isequal(fields{i}, 'wpli')
                ylabel(ax_sub(i),'wpli')
                cmocean('thermal')
            end
            if plotevents
                e = event{theta}{high};
                y = get(gca,'ylim');
                y = y(2);
                x = data.t
                y = y * ones(size(x))-2;
                x(~e) = nan; y(~e) = nan;
                plot(ax_sub(i), x, y, 'Color','k');
            end
            pfc_delta = smoothdata(norm_f(raw{PFC}{delta}));
            %plot(ax_sub(i), data.t, pfc_delta, ...
            %    'Color',[1,1,1], 'linestyle', '--', 'linewidth', 1);
            ylim(ax_sub(i), [0, max(data.f)]);
            %plot(ax_sub, data.t, velocity, ...
            %    'Color',[255,69,0]/255, 'linestyle', ':');
        end
        linkaxes(ax_sub, 'xy');
        xlim([min(data.t), max(data.t)])
        nestable off;

        %% SECTION 2: PLOT TIME!
        % ---------------------------------------------
        % (1) Determine periods of high and low delta
        ax = subplot(2,2,3:4);
        cla
        plotevents = true;
        title(['Session' sprintf('-%s-%d-%d',animal{1}, dayepoch)])
        F = data.f > plot_range(1) & data.f < plot_range(2);
        data.f = data.f(F);
        norm_f  = @(x) ((x / range(quantile(x(:),[0, 0.8])))-0.25)*range(data.f)/2 + min(data.f);
        for field = fields
            data.(field{1}) = data.(field{1})(:,F);
        end
        nestable on;
        ax_sub = gobjects(3,1);
        for i = 1:3
            spectral_data = real(data.(fields{i}));
            ax_sub(i) =  nestplot(3,1,i);
            fields{i}
            if fields{i}(1) == 'w'
                imagesc(ax_sub(i), data.t, data.f, spectral_data', [-0.2, 0.5]);
            elseif fields{i}(1) == 'S'
                %spectral_data = log10(spectral_data+10);
                spectral_data = log10(spectral_data)
                imagesc(ax_sub(i), data.t, data.f, spectral_data', [0, inf]);
                %set(gca,'clim',[-1,2])
            end
            colorbar
            set(gca,'ydir','normal')
            hold on;
            if isequal(fields{i}, 'S1')
                ylabel(ax_sub(i),'CA1')
            elseif isequal(fields{i}, 'S2')
                ylabel(ax_sub(i),'PFC')
            elseif isequal(fields{i}, 'wpli')
                ylabel(ax_sub(i),'wpli')
                cmocean('thermal')
            end
            if plotevents
                e = event{theta}{high};
                y = get(gca,'ylim');
                y = y(2);
                x = data.t
                y = y * ones(size(x))-2;
                x(~e) = nan; y(~e) = nan;
                plot(ax_sub(i), x, y, 'Color','k');
            end
            pfc_delta = smoothdata(norm_f(raw{PFC}{delta}));
            %plot(ax_sub(i), data.t, pfc_delta, ...
            %    'Color',[1,1,1], 'linestyle', '--', 'linewidth', 1);
            ylim(ax_sub(i), [0, max(data.f)]);
            %plot(ax_sub, data.t, velocity, ...
            %    'Color',[255,69,0]/255, 'linestyle', ':');
        end
        linkaxes(ax_sub, 'xy');
        xlim([min(data.t), max(data.t)])
        nestable off;

        % (2) Raw levels
        ax = subplot(3,2,3);
        cla;
        overall_container = nestable('on');

            container = nestplot(1,2,1, 'container', overall_container)
            sca(container)
            title(container, 'Raw Values')
            nestable('on', container)
                %CA1
                nestplot(2,1,1)
                cla
                hold on;
                raw_vals = raw{CA1}{theta}(event{theta}{low})
                h1= histogram(raw_vals(:), 100, 'FaceColor', 'b', 'FaceAlpha', 0.25)
                raw_vals = raw{CA1}{theta}(event{theta}{high})
                h2= histogram(raw_vals(:), 100, 'FaceColor', 'r', 'FaceAlpha', 0.25)
                legend([h1,h2], 'Low','High')
                ylabel('CA1')
                %CA1PFC
                nestplot(2,1,2)
                cla
                hold on;
                raw_vals = raw{CA1PFC}{theta}(event{theta}{low})
                h1= histogram(raw_vals(:), 100, 'FaceColor', 'b', 'FaceAlpha', 0.25)
                raw_vals = raw{CA1PFC}{theta}(event{theta}{high})
                h2= histogram(raw_vals(:), 100, 'FaceColor', 'r', 'FaceAlpha', 0.25)
                legend([h1,h2], 'Low','High')
                ylabel('CA1PFC wpli')
                %CA1PFC
            nestable('off',container)

            container = nestplot(1,2,2, 'container', overall_container)
            sca(container)
            title(container, 'Bootstrap Means')
            nestplot('on', container)
                %CA1
                nestplot(2,1,1, 'container', container)
                cla
                hold on;
                raw_vals = average{CA1}{theta}{low};
                h1= histogram(raw_vals(:), 100, 'FaceColor', 'b', 'FaceAlpha', 0.25)
                raw_vals = average{CA1}{theta}{high};
                h2= histogram(raw_vals(:), 100, 'FaceColor', 'r', 'FaceAlpha', 0.25)
                legend([h1,h2], 'Low','High')
                ylabel('CA1')
                %CA1PFC
                nestplot(2,1,2, 'container', container)
                cla
                hold on;
                raw_vals = average{CA1PFC}{theta}{low};
                h1= histogram(raw_vals(:), 100, 'FaceColor', 'b', 'FaceAlpha', 0.25)
                raw_vals = average{CA1PFC}{theta}{high};
                h2= histogram(raw_vals(:), 100, 'FaceColor', 'r', 'FaceAlpha', 0.25)
                legend([h1,h2], 'Low','High')
                ylabel('CA1PFC')
            nestable off
        nestable off


        % (3) Bootstrap
        ax = subplot(3,2,4);

        % (4) Trigger
        % -----------
        ax = subplot(3,2,5);
        title('Peak triggered windows')
        cla
        nestable on
        x = linspace(-win(1), win(2), winsize);
        ax_sub = gobjects(2,1);
        ax_sub(1) = nestplot(2,1,1)
        avg_low = mean(windows{CA1}{low}, 1);
        hold on
        avg_high = mean(windows{CA1}{high}, 1);
        p1= plot(x, avg_low, 'color', 'b')
        p2=plot(x, avg_high, 'color', 'r')
        h = hline(mean(transform(raw{CA1}{theta})));
        h.Color = 'k'; h.LineStyle = ':';
        ylabel('CA1 Theta')
        %legend([p1,p2], 'low','high')
        ax_sub(2) = nestplot(2,1,2)
        avg_low = mean(windows{CA1PFC}{low}, 1);
        %legend([p1,p2], 'low','high')
        avg_high = mean(windows{CA1PFC}{high}, 1);
        p1= plot(x, avg_low, 'color', 'b')
        hold on
        p2=plot(x, avg_high, 'color', 'r')
        h = hline(mean(raw{CA1PFC}{theta}))
        h.Color = 'k'; h.LineStyle = ':';
        ylabel('ca1-pfc wpli')
        %legend([p1,p2], 'low','high')
        linkaxes(ax_sub, 'x')
        nestable off

        % (5)
        ax = subplot(3,2,6);
        cla
        nestable on
        % (5a) Spectrum in windows
        nestplot(2,2,1)
        difflocs_high = cell(1,nWindows);
        difflocs_low = cell(1,nWindows);
        for window = 1:nWindows
            [peaks, locs] = findpeaks(windows{CA1}{high}(window,:));
            difflocs_low{window} = diff(locs);
            [peaks, locs] = findpeaks(windows{CA1}{low}(window,:));
            difflocs_high{window} = diff(locs);
        end
        difflocs_low = cat(2, difflocs_low{:}) * samprate;
        difflocs_high = cat(2, difflocs_high{:}) * samprate;
        histogram(difflocs_low, 'FaceColor','b','FaceAlpha',0.25); hold on
        histogram(difflocs_high, 'FaceColor','r','FaceAlpha',0.25)
        v =vline( mode(difflocs_high)); v.Color='k'; v.LineStyle = '-';
        set(gca,'yscale','log')
        nestplot(2,2,3)
        difflocs_high = cell(1,nWindows);
        difflocs_low = cell(1,nWindows);
        for window = 1:nWindows
            [peaks, locs] = findpeaks(windows{CA1PFC}{high}(window,:));
            difflocs_low{window} = diff(locs);
            [peaks, locs] = findpeaks(windows{CA1PFC}{low}(window,:));
            difflocs_high{window} = diff(locs);
        end
        difflocs_low = cat(2, difflocs_low{:}) * samprate;
        difflocs_high = cat(2, difflocs_high{:}) * samprate;
        histogram(difflocs_low, 'FaceColor','b','FaceAlpha',0.25); hold on
        histogram(difflocs_high, 'FaceColor','r','FaceAlpha',0.25)
        v =vline( mode(difflocs_high)); v.Color='k'; v.LineStyle = '-';
        % (5b) Peak-to-peak histogram in windows
        nestplot(2,2,2);
        nestplot(2,2,4)
        nestable off
    end
end
