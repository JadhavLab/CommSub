% RY: July 18 made a change to resolve nan issue in some animasl

clear results
resolution = 50;
strongHeadDirection = true;
ProgressBar.deleteAllTimers()
acnt = 0;
animals = {'ER1','KL8','JS12','JS13','JS14','JS15','JS17','JS21','ZT2'};
for animal = progress(animals,'Title','animal')
    acnt=acnt+1;

    animalToPath(animal{1});
    load( string(animal) + "linpos01");

    for epoch = progress(2:2:16,'Title','epoch')

        %try
        for traj = 1:4

            %% Lindist
            lp = linpos{1}{epoch};
            load( string(animal) + "cgramcnew" + sprintf('-01-%02d.mat', epoch));
            c = ffend(cgramcnew);

            trajinds = lp.statematrix.traj == traj;
            headinds = abs(max(lp.statematrix.segmentHeadDirection,[],2))>0.5;
            if strongHeadDirection
                inds = headinds & trajinds;
            else
                inds = trajinds;
            end
            inds = find(inds);
                
            timeTranslate = lookup(lp.statematrix.time(inds), c.t);
            c.t = c.t(timeTranslate);
            c.C = c.C(timeTranslate,:);
            c.S1 = c.S1(timeTranslate,:);
            c.S2 = c.S2(timeTranslate,:);
            c.wpli = c.wpli(timeTranslate,:);

            bins = discretize(lp.statematrix.lindist,resolution);
            bins = bins(inds);

            binmapping = findgroups(bins);
            B = [unique(binmapping), unique(bins)];

            results(acnt).animal = string(animal);
            for field = {'C','S1','S2','wpli'}
                Y = splitapply(@(x) nanmean(x,1), c.(field{1}), binmapping);
                Ynew = nan(resolution, numel(c.f));
                Ynew(B(:,2),:) = Y(:,:);
                %Y = smoothdata(Y,2);
                %Y = smoothdata(Y,1);
                results(acnt, epoch/2, traj).(field{1}) = Ynew;
            end
            results(acnt, epoch/2, traj).t = c.t;
            results(acnt, epoch/2, traj).f = c.f;
            
            doWellTimes = false;
            if doWellTimes
                %% Well times
                
                % Compute time to nearest well event
                time = lp.statematrix.lindist;
                welltime = interp1(time_of_wellevent, time_of_well_event, time, 'nearest'); % map each time to nearest well event
                welltime=time-welltime;
                % Discretize!
                edges = linspace(-3,3,50);
                bins = discretize(welltimes, edges);

                binmapping = findgroups(bins);
                B = [unique(binmapping), unique(bins)];

                resultsT(acnt).animal = string(animal);
                for field = {'C','S1','S2','wpli'}
                    Y = splitapply(@(x) nanmean(x,1), c.(field{1}), binmapping);
                    Ynew = nan(resolution, numel(c.f));
                    Ynew(B(:,2),:) = Y(:,:);
                    %Y = smoothdata(Y,2);
                    %Y = smoothdata(Y,1);
                    resultsT(acnt, epoch/2, traj).(field{1}) = Ynew;
                end
                resultsT(acnt, epoch/2, traj).t = c.t;
                resultsT(acnt, epoch/2, traj).f = c.f;
            end
            
        end

        %catch ME
        %    disp(ME)
        %end

    end

end
