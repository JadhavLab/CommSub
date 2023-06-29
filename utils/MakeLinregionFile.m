%% Create Linearized Region File
% Which animals to use?

%animal={'JS17','JS21','JS15','HPa','HPb','HPc','YD6','ER1','KL8','JS12','JS13','JS14'};
animal={'ER1'};
nregion = 10;
%% 
% Now we iterate over them and create the linearized region file which is 
% like trajinfo but includes information about the time by time linear region.

for A = animal
a = animaldef(A{1});

switch A{1}
    case {'HPa','HPb','HPc'}, days=1:8;
    case {'YD6','ER1','KL8','JS12','JS13','JS14','JS15','JS17','JS21'}, days=1;
end

linpos=loaddatastruct(a{2:3},'linpos',days);
pos=loaddatastruct(a{2:3},'pos',days);
[fetch] = cellfetch(pos,'data')
vel = DFTFry_getvelpos(a{2:3},fetch.index);
trajinfo=loaddatastruct(a{2:3},'perftrajinfo',days);
try 
    dio=loaddatastruct(a{2:3},'DIO',days);
catch
    dio=loaddatastruct(a{2:3},'dio',days);
end

% We will populate a number of informations per each trajectory.
% ----- ----- ----- ----- -----
% (1) Region, Traj
% (2) Regional ExclusionTime
% (3) DIO Corrected intertrial boundary

clear out
nDays = numel(trajinfo);
for d = 1:nDays
nEps = numel(trajinfo{d});
for e = 1:nEps
    fprintf('Day %d Epoch %d\n',d,e);
    if isempty(linpos{d}{e});continue;end
    
    [times,normlindist]=relevent_linpos_features(linpos{d}{e});
    lregion{d}{e}.time = times;
    lregion{d}{e}.desc = sprintf("TR->trajectory region; ");
    lregion{d}{e}.TR = [];
    lregion{d}{e}.TB = [];
    lregion{d}{e}.tinctimes=[];
    
    if isempty(trajinfo{d}{e});continue;end
    % Compute (1)
    for trajs = 1:size(trajinfo{d}{e}.trajtime,1)
        ttime = times >= trajinfo{d}{e}.trajtime(trajs,1) & times < trajinfo{d}{e}.trajtime(trajs,2);
        if trajs == size(trajinfo{d}{e}.trajtime,1); ttime(end)=true; end
        lregion{d}{e}.TR=[lregion{d}{e}.TR; repmat(trajs,sum(ttime),1), ceil(normlindist(ttime)*nregion)];
        % (3) Exclusion times for trajs (these can be combined with the rexctimes via combineExlcusionPeriods)
        lregion{d}{e}.tinctimes(trajs,:) = trajinfo{d}{e}.trajtime(trajs,:);
        % REWARDED AND TRAJBOUND
        lregion{d}{e}.TB =[lregion{d}{e}.TB; ...
            repmat(trajinfo{d}{e}.rewarded(trajs),sum(ttime),1), repmat(trajinfo{d}{e}.trajbound(trajs),sum(ttime),1),...
            ];
        % PERFORMANCE
        lregion{d}{e}.tperf = trajinfo{d}{e}.pcorrect; % trial-by-trial pcorrect estimates
        lregion{d}{e}.epochperf = [trajinfo{d}{e}.avgoutbound; trajinfo{d}{e}.avginbound]; % epoch pcorrect estimates (mean through the session)
    end
    
    
    lregion{d}{e}.TR(lregion{d}{e}.TR(:,2)==0,2) = 1;
    lregion{d}{e}.trajbound = trajinfo{d}{e}.trajbound;
    lregion{d}{e}.rewarded = trajinfo{d}{e}.rewarded;
    
    lregion{d}{e}.XY  = [pos{d}{e}.data(:,2:3)];
    lregion{d}{e}.vel = [vel{d}{e}.absvel(:)];

    
    % Correcting for rare case in which linpos has 2-3 more time
    % time points than pos
    plot_check  = false;
    if size(pos{d}{e}.data,1) ~= numel(times)
        pos_time = vel{d}{e}.time;
        
        if plot_check
            figure(1);clf;
            subplot(121)
            plot(repmat(pos_time,[1 2]), lregion{d}{e}.XY);hold on;
            subplot(122)
            plot(pos_time, lregion{d}{e}.vel);hold on;
        end
        
        lregion{d}{e}.XY = interp1(pos_time,lregion{d}{e}.XY,times,'cubic','extrap');
        lregion{d}{e}.vel = interp1(pos_time,lregion{d}{e}.vel,times,'cubic','extrap');
        
        if plot_check
            figure(1);
            ax1=subplot(121)
            plot(repmat(times,[1 2]), lregion{d}{e}.XY);
            axis off;
            ax2=subplot(122)
            plot(times, lregion{d}{e}.vel);
            axis off;
            linkaxes([ax1,ax2],'x')
            set(findobj('type','line'),'linewidth',2)
        end
        
    end
    assert( numel(lregion{d}{e}.vel) == numel(lregion{d}{e}.time) );
    lregion{d}{e}.lindist = normlindist;
    if isfield(linpos{d}{e}.statematrix,'XY_proj')
        lregion{d}{e}.XY = [lregion{d}{e}.XY , linpos{d}{e}.statematrix.XY_proj];
    end
    
    % (2) Compute exclusion times per region
    lregion{d}{e}.rinctimes=cell(1,nregion);
    for r = 1:nregion
        lregion{d}{e}.rinctimes{r} = getIncludePeriods(times,lregion{d}{e}.TR(:,2)==r);
    end
    
    % Make tperf timewise
    clear trial_of_time
    tr_pick = bsxfun(@lt,lregion{d}{e}.time,lregion{d}{e}.tinctimes(:,2)'); % first, extract all possible trial start times that are less than the spike times
    for t = 1:size(tr_pick,1)
        trial = find(tr_pick(t,:),1,'first');
        if ~isempty(trial)
            trial_of_time(t) = trial; % trial is the last where the inequality holds true
        else
            trial_of_time(t) = nan;
        end
    end
    nans = find(isnan(trial_of_time));
    trial_of_time(nans) = trial_of_time(max(nans-1,1));
    lregion{d}{e}.tperf_timewise = lregion{d}{e}.tperf(trial_of_time,1);
    
    lregion{d}{e}.date = date;
    
end
end

savedatastruct(lregion,a{2:3},'lregion');

end

%% Helper Functions
%%
function [times,normlindist] = relevent_linpos_features(linpos)
    times = linpos.statematrix.time;
    if isfield(linpos.statematrix,'normlindist')
        normlindist=linpos.statematrix.normlindist;
    else
        normlindist=linpos.statematrix.lindist/max(linpos.statematrix.lindist);
    end
end
