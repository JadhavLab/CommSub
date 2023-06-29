% RY: July 18 made a change to resolve nan issue in some animasl
clear results
resolution = 50;
strongHeadDirection = true;
binField = "normPathLength";
ProgressBar.deleteAllTimers()
acnt = 0;
animals = ["RY16"]
for animal = progress(animals, 'Title', 'animal')

    acnt=acnt+1;

    for day   = progress(36,'Title','epoch')
    for epoch = progress([2 4 5 7],'Title','epoch')

        %% Load
        E = ndb.load(animal, "egocentric", 'indices',  [day epoch], 'get', true);
        C  = ndb.load(animal, "cgramcnew",...
            'indices', [day  epoch], 'get', true);
        E.normPathLength = E.currentPathLength./E.currentPathLengthMax;

        [dups] = util.getduplicates(C.t);
        C.t(dups) = [];
        C.C(dups,:) = [];
        C.S1(dups,:) = [];
        C.S2(dups,:) = [];
        C.wpli(dups,:) = [];

        for cuemem  = 0:1
        for region  = 0:1
        for correct = 0:1

            %  Compute all times matching our properties
            relInds = ~isnan(E.(binField)) & E.cuemem == cuemem &...
                      E.correct == correct &...
                      E.region == region;
            relInds = find(relInds);
            relTimes = E.postime(relInds);

            % Remove any duplicate times
            c = C;
                
            % Grab all times that match our properties
            timeTranslate = interp1(c.t, 1:numel(c.t), relTimes, 'nearest');
            timeTranslate = timeTranslate(~isnan(timeTranslate));
            c.t    = c.t(timeTranslate);
            c.C    = c.C(timeTranslate,:);
            c.S1   = c.S1(timeTranslate,:);
            c.S2   = c.S2(timeTranslate,:);
            c.wpli = c.wpli(timeTranslate,:);

            % Bin the relevent times
            discProperty  = E.(binField);
            bins = discretize(discProperty, resolution);
            bins = bins(relInds);
            if isempty(bins)
                continue
            end

            binmapping = findgroups(bins);
            B = [unique(binmapping), unique(bins)];

            results(acnt).animal = string(animal);
            for field = {'C','S1','S2','wpli'}
                Y = splitapply(@(x) nanmean(x,1), c.(field{1}), binmapping);
                Ynew = nan(resolution, numel(c.f));
                Ynew(B(:,2),:) = Y(:,:);
                results(acnt, epoch, cuemem+1, region+1, correct+1).(field{1}) = Ynew;
            end
            results(acnt, epoch, cuemem+1, region+1, correct+1).t = c.t;
            results(acnt, epoch, cuemem+1, region+1, correct+1).f = c.f;
            
            % Get the start and end of our property of interest and let's
            % grab -3 seconds before and +3 seconds after
            doStartStop = false;
            if doStartStop
                %% Compute padded nan-edge times
                
                % Compute distance of each padded nan-edge time to each start-well time
                time = paddedNanEdgeTimes;
                N = [0; diff(isnan(discProperty))];
                S(:,1) = find(N);
                S(:,2) = discProperty(max(1, find(N) + 2)) > 0.1;
                S(:,3) = discProperty(min(1, find(N) - 2)) < 0.1;
                starttimes = time(S(S(:,2) == 1,1));
                stoptimes  = time(S(S(:,3) == 1,1));
                startwelltime = interp1(starttimes, starttimes, time, 'nearest'); % map each time to nearest well event
                stopwelltime = interp1(starttimes,  starttimes, time, 'nearest'); % map each time to nearest well event
                startwelltime = time - startwelltime;
                stopwelltime = time - stopwelltime;

                % Set each 

                % Discretize!
                startedges = linspace(-3, 0.25, 25);
                startbins  = discretize(startwelltime, startedges);
                stopedges = linspace(-0.25, 3, 25);
                stopbins  = discretize(stopwelltime, stopedges);

                % START 
                binmapping = findgroups(startbins);
                B = [unique(binmapping), unique(startbins)];
                if any(isnan(binmapping))
                    dummyBin = max(binmapping) + 1;
                    binmapping(isnan(binmapping)) = dummyBin;
                end

                resultsStart(acnt).animal = string(animal);
                for field = {'C','S1','S2','wpli'}
                    Y              = splitapply(@(x) nanmean(x,1), c.(field{1}), binmapping);
                    if ~isnan(dummyBin)
                        Y(binmapping == dummyBin,:) = [];
                    end
                    Ynew           = nan(resolution, numel(c.f));
                    Ynew(B(:,2),:) = Y(:,:);
                    %Y = smoothdata(Y,2);
                    %Y = smoothdata(Y,1);
                    resultsStart(acnt, epoch, cuemem+1, region+1, correct+1).(field{1}) = Ynew;
                end
                resultsT(acnt, epoch, cuemem+1, region+1, correct+1).t = c.t;
                resultsT(acnt, epoch, cuemem+1, region+1, correct+1).f = c.f;

                % STOP
                binmapping = findgroups(stopbins);
                B = [unique(binmapping), unique(stopbins)];

                resultsStop(acnt).animal = string(animal);
                for field = {'C','S1','S2','wpli'}
                    Y              = splitapply(@(x) nanmean(x,1), c.(field{1}), binmapping);
                    Ynew           = nan(resolution, numel(c.f));
                    Ynew(B(:,2),:) = Y(:,:);
                    %Y = smoothdata(Y,2);
                    %Y = smoothdata(Y,1);
                    resultsStop(acnt, epoch/2, traj).(field{1}) = Ynew;
                end
                resultsStop(acnt, epoch/2, traj).t = c.t;
                resultsStop(acnt, epoch/2, traj).f = c.f;
            end
            
        end

        %catch ME
        %    disp(ME)
        %end

    end
    end
    end

end
