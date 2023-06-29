% RY: July 18 made a change to resolve nan issue in some animasl

clear results
resolution = 50;
fieldset = {'C','S1','S2','wpli','plv','phi'};
strongHeadDirection = true;
% traj_dims_I_care_about = ["trajbound", "rewarded", "choice"];
% ProgressBar.deleteAllTimers()

% Determine the animal set
animal_set = "mingxin";
if animal_set == "mingxin"
    animals = {'SM8'};
    epochs = 2:2:4;
    traj_dims_I_care_about = ["trajbound", "choice"];
    maze="T";
    %     task = 'Tmaze';
elseif animal_set == "justin"
    animals = {'ER1','KL8','JS12','JS13','JS14','JS15','JS17','JS21','ZT2'};
    epochs = 2:2:16;
    traj_dims_I_care_about = ["trajbound"];
    maze="W";
end

acnt = 0;
for animal = progress(animals,'Title','animal')
    acnt=acnt+1;
    
    animalToPath(animal{1});
    linpos = ndb.load(animal{1}, 'linpos');
    indices = ndb.indicesMatrixForm(linpos);
    
    trajinfo = ndb.load(animal{1}, 'trajinfo');
    
    for epoch = progress(epochs,'Title','epoch')
        
        %try
        %         for traj = 1:4
        
        %% Lindist
        lp = linpos{1}{epoch};
        ti_epoch = trajinfo{1}{epoch};
        trajtime = trajinfo{1}{epoch}.trajtime;
        time = lp.statematrix.time;
        startinds = lookup(trajtime(:,1), lp.statematrix.time);
        endinds = lookup(trajtime(:,2), lp.statematrix.time) - 1;
        load( string(animal) + "cgramcnew" + sprintf('-01-%02d.mat', epoch));
        c = ffend(cgramcnew);
        
        
        fields_I_care_about = arrayfun(@(field) ti_epoch.(field), traj_dims_I_care_about, 'UniformOutput', false);
        T = table(fields_I_care_about{:}, 'VariableNames', traj_dims_I_care_about);
        
        group_item_values = cell(1, numel(traj_dims_I_care_about));
        [groups, group_item_values{:}] = findgroups(fields_I_care_about{:});
        
        uGroups = unique(groups);
        
        for g = uGroups'
            
            
            % Get all times that match our group
            HIGH = 2; LOW = 1;
            matching_times = false(1,numel(lp.statematrix.time));
            these_trajectories = trajtime(groups == g,:);
            for trajectory = 1:size(these_trajectories, 1)
                matching_times = matching_times(:) | ...
                    (time(:) >= these_trajectories(trajectory, LOW) & time(:) < these_trajectories(trajectory, HIGH));
            end
            
            headinds = abs(max(lp.statematrix.segmentHeadDirection,[],2))>0.5;            
            if strongHeadDirection
                inds = matching_times(:) & headinds(:);
%                 inds = headinds & trajinds;
            else
                inds = matching_times(:);
%                 inds = trajinds;
            end
            
            inds = find(inds);
            
            timeTranslate = lookup(lp.statematrix.time(inds), c.t);
            c.t = c.t(timeTranslate);
            for field = string(fieldset)
                c.(field)    = c.(field)(timeTranslate,:);
            end
            
            bins = discretize(lp.statematrix.lindist,resolution);
            bins = bins(inds);
            
            binmapping = findgroups(bins);
            B = [unique(binmapping), unique(bins)];
            
            % Compute address where to put results
            address = cellfun(@(address_item) address_item(g) + 1, group_item_values);
            address = num2cell(address);
            
            % Compute results
            results(acnt).animal = string(animal);
            for field = fieldset
                if field == "phi"
                    func = @(x) circ_mean(x);
                else
                    func = @(x) nanmean(x,1);
                end
                Y = splitapply(func, c.(field{1}), binmapping);
                Ynew = nan(resolution, numel(c.f));
                Ynew(B(:,2),:) = Y(:,:);
                %Y = smoothdata(Y,2);
                %Y = smoothdata(Y,1);
%                 results(acnt, epoch/2, traj).(field{1}) = Ynew;
                results(acnt, epoch/2, address{:}).(field{1}) = Ynew;
                %if field == "phi"
                %    fig(animal + " epoch_" + epoch +  " group_" + g); imagesc(1:50,c.f,Ynew'); axis xy; crameri bamO; colorbar;  xlabel("Linear Bin"); ylabel("Freq (Hz)"); sgtitle(get(gcf,'Name'));
                %end
            end
%             results(acnt, epoch/2, traj).t = c.t;
%             results(acnt, epoch/2, traj).f = c.f;
            results(acnt, epoch/2, address{:}).t = c.t;
            results(acnt, epoch/2, address{:}).f = c.f;
            
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
                for field = fieldset
                    Y = splitapply(@(x) nanmean(x,1), c.(field{1}), binmapping);
                    Ynew = nan(resolution, numel(c.f));
                    Ynew(B(:,2),:) = Y(:,:);
                    %Y = smoothdata(Y,2);
                    %Y = smoothdata(Y,1);
%                     resultsT(acnt, epoch/2, traj).(field{1}) = Ynew;
                    results(acnt, epoch/2, trial).(field{1}) = Ynew;
                end
%                 resultsT(acnt, epoch/2, traj).t = c.t;
%                 resultsT(acnt, epoch/2, traj).f = c.f;
                results(acnt, epoch/2, trial).t = c.t;
                results(acnt, epoch/2, trial).f = c.f;
            end
        end
    end
    %catch ME
    %    disp(ME)
    %end
    
% end

% Save
results = nd.dimLabel(results, 1:(2+length(traj_dims_I_care_about)), ["animal", "epoch", traj_dims_I_care_about]);
results_ndb = nd.toNdBranch(results);
ndb.save(results_ndb{acnt}, animal{1}, 'linearizedspecnew');

end

