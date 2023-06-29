function globalripple = generateGlobalRipple(animalprefix)

days = 1;
minstd = 3;
minrip = 1;
minenergy = 0;
velfilter = 4; %velthresh
matcheegtime = 0;

%%
tetinfo = strcat(animalprefix,'tetinfo.mat'); % get tetrode info
load(tetinfo);
[eps,~] = getRunningSessions(animalprefix); %only using the sleep sessions.

%%
clear index;
for day = days
    nn = 0;
    index = zeros(length(eps),2);
    
    for ep = eps
        nn = nn+1;
        index(nn,1) = day;
        index(nn,2) = ep;
    end
    
    for id = 1:nn
        
        d = index(id,1);
        e = index(id,2);
        tetfilter = 'isequal($descrip, ''riptet'')';
        tetlist =  evaluatefilter(tetinfo{d}{e}, tetfilter);
        tetlist = unique(tetlist(:,1))';
        load(sprintf('%sripples%02d', animalprefix, d));
        
        % Get the exclude periods
        % ------------------------
        if mod(e,2) == 0 
            % Run epochs: exclude periods are high velocity periods
            if ~isempty(velfilter)
                load(sprintf('%spos%02d.mat',animalprefix, d)); % get animal [time, positions, velocity]
                if size(pos{d}{e}.data,2) > 5  % get velocity
                    velocity = pos{d}{e}.data(:,9);% smoothed velocity
                else
                    velocity = pos{d}{e}.data(:,5);
                end
                postimes = pos{d}{e}.data(:,1);  % get time stamps
                immobile = vec2list(velocity < velfilter, postimes); % generate [start end] list of immobile epochs
                excludeperiods = immobile;
            end
        else
        % Sleep epochs: exclude periods are non-SWS periods
        animaldir = "C:\Users\BrainMaker\commsubspace\SingleDayExpt\"+animalprefix+"_direct\"
             load(sprintf('%s%ssws0%d.mat',animaldir,animalprefix,d)); % get sws time
             swsep = sws{d}{e};
             swslist(:,1) = swsep.starttime;
             swslist(:,2) = swsep.endtime;
             excludeperiods = swslist;
        end
        
        r = ripples{d}{e}{tetlist(1)};
        times = r.timerange(1):0.001:r.timerange(end);
        nrip = zeros(size(times)); nstd=[];
        ripplestd = zeros(size(times));
        for t = 1:length(tetlist)
             tmprip = ripples{d}{e}{tetlist(t)};
             % get the indeces for the ripples with energy above minenergy
             % and maxthresh above minstd
             rvalid = find((tmprip.energy >= minenergy) & (tmprip.maxthresh >= minstd) & (isExcluded(tmprip.midtime,excludeperiods)));
             rtimes = [tmprip.starttime(rvalid) tmprip.endtime(rvalid)];
             tmpripplestd = [tmprip.maxthresh(rvalid) tmprip.maxthresh(rvalid)];
             tmprippleenergy = [tmprip.energy(rvalid) tmprip.energy(rvalid)];
             % create another parallel vector with bordering times for zeros
             nrtimes = [(rtimes(:,1) - 0.00001) (rtimes(:,2) + 0.00001)];
             rtimes = reshape(rtimes', length(rtimes(:)), 1);
             rtimes(:,2) = 1;
             tmpriplestd = [rtimes(:,1) tmpripplestd(:)];
             nrtimes = [r.timerange(1) ; reshape(nrtimes', ...
                 length(nrtimes(:)), 1) ; r.timerange(2)];
             nrtimes(:,2) = 0;
             % create a new list with all of the times in it
             tlist = sortrows([rtimes ; nrtimes]);
             [junk, ind] = unique(tlist(:,1));
             tlist = tlist(ind,:);

             stdlist = sortrows([tmpriplestd ; nrtimes]);
             stdlist =stdlist(ind,:);
             nrip = nrip + interp1(tlist(:,1), tlist(:,2), times, 'nearest');
             nstd(t,:) = interp1(stdlist(:,1), stdlist(:,2), times, 'nearest');  % carry forward amplitude of ripple
        end
        %find the start and end borders of each ripple
        inripple = (nrip >= minrip);
        startrippleind = find(diff(inripple) == 1)+1;
        endrippleind = find(diff(inripple) == -1)+1;
        ripplestdout = [];
        if ((length(startrippleind) > 1) && (length(endrippleind) > 1))
            if (endrippleind(1) < startrippleind(1))
                endrippleind = endrippleind(2:end);
            end
            if (endrippleind(end) < startrippleind(end))
                startrippleind = startrippleind(1:end-1);
            end
            startripple = times(startrippleind);
            endripple = times(endrippleind);

            % Get amplitude of "global" ripple: maximum across tetrodes
            max_nstd = max(nstd,[],1);
            ampripple = max_nstd(startrippleind);
    %         for i = 1:length(startripple)
    %             ripplestdout(i,1) = max(ripplestd(startrippleind(i):endrippleind(i)));
    %         end

            out = [startripple(:) endripple(:) ampripple(:)]; % amplitude of ripple

            out = out( ((out(:,2)-out(:,1))>.050),:);  % ripples separated by 50 ms
        else
            out = [];
            ripplestdout = [];
        end
        globalripple{d}{e} = out;
        clear out
        clear swslist
    end
%     save(sprintf('data/JS15/%sglobalripple%02d.mat', animalprefix, d), 'globalripple');
%     clear globalripple
end

end

