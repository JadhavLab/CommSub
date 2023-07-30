function [behavior, times_to_throwout, unique_times_inds] = ...
    lookup(animal, times, inds, varargin)
% [behavior, times_to_throwout, unique_times] = ...
%    lookup(animal, times, inds, varargin)
%
% Given a set of times, looks up everything relevant about the behavior and
% outputs into table form for introspection with gramm
%
% Inputs
% ------
% animal : string
%     Animal name
% times : array
%     Nx1 array of times to look up
% inds : array
%     Nx2 array of indices to look up
% varargin : key/value pairs
%     'valueOnly' : bool
%         Whether to not return ccolumns that function as indices (time, traj, day, epoch)
%     'scaleVars' : bool
%         Whether to scale vars
%     'throwOutSleep' : bool
%         Whether to throw out sleep
%     'throwOutThresh' : int
%         If 2 seconds from nearest point, toss it
%
% Outputs
% -------
% behavior : table
%     Table of behavior
% times_to_throwout : array
%     Times that were thrown out
% unique_times : array
%     Unique times that were looked up

disp("Starting lookup for " + animal);
tic;

ip = inputParser;
ip.addParameter('valueOnly',false); % Whether to not return ccolumns that function as indices (time, traj, day, epoch)
ip.addParameter('scaleVars',false); % Whether to scale vars
ip.addParameter('throwOutSleep',true); % Whether to scale vars
ip.addParameter('throwOutThresh',2); % If 2 seconds from nearest point, toss it
ip.parse(varargin{:})
Opt = ip.Results;

if nargin < 2
    times = [];
end
if nargin < 3
    inds = [];
end

if Opt.throwOutSleep
    % Section being developed
    task    = ndb.load(animal,'task');
    run_inds = evaluatefilter(task, 'string($type) == "run"');
    if ~isempty(inds)
        inds = intersect(inds, run_inds, 'rows');
    else
        inds = run_inds;
    end
end

pos     = ndb.load(animal, 'pos',     'indices', inds);
linpos  = ndb.load(animal, 'linpos',  'indices', inds);
lregion = ndb.load(animal, 'lregion', 'indices', inds);
dayep   = ndb.indicesMatrixForm(lregion);
lregion = ndb.toNd(lregion);
pos     = ndb.toNd(pos);
linpos  = ndb.toNd(linpos);

% Crate label fields;
% Drop all irrelevent fields
irrelevent = ["tinctimes", "tperf", "epochperf", "trajbound", "rewarded", "rinctimes", "date", "desc"];
lregion = rmfield(lregion, irrelevent);

% Rename and munge fields again
% -----------------------------
for field = ["region","traj","trajbound","rewarded","X","Y"]
    lregion = arrayfun(@(lregion) setfield(lregion, field, []), lregion);
end
lregion = nd.apply(lregion, @(x) setfield(x, 'region',    x.TR(:, 2)));
lregion = nd.apply(lregion, @(x) setfield(x, 'traj',      x.TR(:, 1)));
lregion = nd.apply(lregion, @(x) setfield(x, 'trajbound', x.TB(:, 2)));
lregion = nd.apply(lregion, @(x) setfield(x, 'rewarded',  x.TB(:, 1)));
lregion = nd.apply(lregion, @(x) setfield(x, 'X',         x.XY(:, 1)));
lregion = nd.apply(lregion, @(x) setfield(x, 'Y',         x.XY(:, 2)));
irrelevent = ["TB", "TR","XY"];
lregion = rmfield(lregion, irrelevent);
lregion = nd.dimLabel(lregion, 1:2, ["day","epoch"]);
try
    pos = arrayfun(@(pos) setfield(pos, 'headdir', []), pos);
    pos = nd.apply(pos, @(x) setfield(x, 'headdir',    x.data(:, 4)));
    headdir = nd.fieldGetCell(pos, 'headdir');
    headdir = cat(1, headdir{:});
catch E
end

% Concatonate into table
lregion = tidyData.fromNd(lregion);
lregion = lregion((lregion.day==1),:); % HARDCODING DAY TO 1
% add headdir
try
    lregion.headdir = headdir((lregion.day==1),:);
catch E
end
% add trajclass (1-left out, 2-left in, 3 right out, 4 right in)
linpos    = nd.unnest(linpos, 'statematrix', 'traj');
trajclass = nd.fieldGetCell(linpos, 'traj');
trajclass = cat(1, trajclass{:});
%try
% 
lregion.trajclass = trajclass;
lregion.leftright = ismember(trajclass, [1, 2]); % literally is it a left traj
%  Compute the upcoming and previous trajectory class
lregion.future    = repmat(nan, height(lregion), 1);
lregion.previous  = repmat(nan, height(lregion), 1);
lregion.futureRewarded    = repmat(nan, height(lregion), 1);
lregion.previousRewarded  = repmat(nan, height(lregion), 1);
[groups, epochs, trajs] = findgroups(categorical(lregion.epoch), categorical(lregion.traj));
uGroups = unique(groups);
for group = uGroups(:)'

    epoch = int16(epochs(group));
    traj  = int16(trajs(group));
    pastTraj   = traj-1;
    futureTraj = traj+1;
    past   = lregion.traj == pastTraj | lregion.epoch == epoch;
    future = lregion.traj == futureTraj | lregion.epoch == epoch;
    fillout = @(x, y) repmat(median(y), sum(x), 1);

    if ~isempty(past)
        lregion.past(groups==group) = fillout(groups==group, lregion.trajclass(past));
        lregion.pastRewarded(groups==group) = fillout(groups==group, lregion.rewarded(past));
        lregion.pastLeftRight(groups==group) = fillout(groups==group, lregion.leftright(past));
    end
    if ~isempty(future)
        lregion.future(groups==group) = fillout(groups==group, lregion.trajclass(future));
        lregion.futureRewarded(groups==group) = fillout(groups==group, lregion.rewarded(future));
        lregion.futureLeftRight(groups==group) = fillout(groups==group, lregion.leftright(future));
    end
end
%catch E
    %keyboard
%end
% Add individual in/outbound
tperf_inbound  = lregion.tperf_timewise(lregion.trajbound==1);
time_inbound  = lregion.time(lregion.trajbound==1);
tperf_outbound = lregion.tperf_timewise(lregion.trajbound==0);
time_outbound  = lregion.time(lregion.trajbound==0);
lregion.tperf_inbound = interp1(time_inbound, tperf_inbound, lregion.time, 'previous'); % could be a goood case for next or previous here
lregion.tperf_outbound = interp1(time_outbound, tperf_outbound, lregion.time, 'previous'); % could be a goood case for next or previous here

% Add derivative of in/outbound and all
for field = ["tperf_inbound","tperf_outbound"]
    der.(field) = [0; diff(lregion.(field))];
    der.(field)(der.(field) == 0) = nan;
    der.(field) = fillmissing(der.(field), 'next');
    lregion.(replace(replace(field,'tperf','deriv'),'timewise_','')) = der.(field);
end
lregion.deriv = zeros(height(lregion),1);
lregion.deriv(lregion.trajbound == 0) = lregion.deriv_outbound(lregion.trajbound == 0);
lregion.deriv(lregion.trajbound == 1) = lregion.deriv_inbound(lregion.trajbound == 1);
lregion.accel = smoothdata([0; diff(lregion.vel)]);

lregion.directional_lindist = pi * lregion.lindist;
lregion.directional_lindist(lregion.trajbound == 1) = lregion.directional_lindist(lregion.trajbound == 1) + (pi - lregion.directional_lindist(lregion.trajbound == 1));
clear i
lregion.directional_lindist_imag = exp(i * lregion.directional_lindist);

% Calculate measures of behavioral variablilty : zIdPhi and moving-traj-win velocity variance
% -------------------------------------------------------------------------------------------
resolution     = 40;
grid_of_ldist  = linspace(0, 1, resolution + 1);
lregion.ldistB = discretize(lregion.lindist, grid_of_ldist);

try
    % -------------
    lregion.trajall       = findgroups(lregion.epoch, lregion.traj);
    [groups, traj, ldist] = findgroups(lregion.trajall, lregion.ldistB);
    uGroups = unique(groups);
    idphi = zeros(numel(traj), numel(ldist));
    for g = uGroups'

        headdir = lregion(g,:).headdir;
        vel = lregion(g,:).vel;
        dphi = abs(diff(headdir));
        idphi(traj(g), ldist(g)) = sum(dphi);
        dV = abs(diff(vel));
        idv(traj(g), ldist(g)) = sum(dV);

    end
catch E
end
% Lookup times
% 
if isempty(times)
    inds = true(size(lregion.time));
else
    inds = interp1(lregion.time, 1:height(lregion), times, 'nearest');
end
times_to_throwout = isnan(inds);

if Opt.throwOutThresh > 0
    disp("Throwing out times greater than " + Opt.throwOutThresh + " seconds away");
    chunksize = 5000;
    for t = progress(1:chunksize:numel(times),'Title','Throwing out bad times')
        stop = min(t+chunksize-1,numel(times));
        tinds = t:stop;
        time = times(tinds);
        toss_inds = all( abs( lregion.time - time(:)' ) > Opt.throwOutThresh,  1);
        times_to_throwout(tinds(toss_inds)) = true;
    end
end

% Return behavior
inds(times_to_throwout) = [];
behavior = lregion(inds,:);

% OPtional transformations
if Opt.valueOnly
    behavior.time   = [];
    behavior.day    = [];
    behavior.epoch  = [];
    behavior.traj   = [];
    behavior.region = [];
    behavior.X      = [];
    behavior.Y      = [];
end
if Opt.scaleVars
    for field = string(fieldnames(behavior))'
        try
            behavior.(field) = behavior.(field) ./ sqrt(nanvar(behavior.(field)));
        catch
        end
    end
end

% ------------------------------------
%% Ziyi's addToBehavior function items
% ------------------------------------
[~, unique_times_inds, ~] = unique(behavior.time);
behavior = behavior( unique_times_inds,:);

%% identify the decision, error, reward times to the table
behavior.wellTimes = (behavior.trajbound == 0  & behavior.lindist > 0.99) | ...
                     (behavior.trajbound == 1 & behavior.lindist <0.01) ; %0.01 inbound
behavior.rewardTimes = behavior.wellTimes & behavior.rewarded ==1;
behavior.outBounderrorTimes = behavior.wellTimes & behavior.rewarded ==0;
behavior.inBounderrorTimes = behavior.trajbound == 1 & behavior.lindist > 0.99  & behavior.rewarded ==0;
behavior.errorTimes = behavior.outBounderrorTimes | behavior.inBounderrorTimes;

behavior.outBoundChoiceTimes = behavior.trajbound == 0 & behavior.lindist >= 0.2 & behavior.lindist <= 0.4;
behavior.inBoundChoiceTimes = behavior.trajbound == 1 & behavior.lindist >= 0.5 & behavior.lindist <= 0.6;

%% interp the inbound and outbound performances to all times
tperfInbound = (behavior.tperf_timewise(behavior.trajbound ==1));
tperfOutbound = (behavior.tperf_timewise(behavior.trajbound ==0));

behavior.tperfInbound = interp1(behavior.time(behavior.trajbound ==1),...
               tperfInbound, behavior.time);
           
behavior.tperfOutbound = interp1(behavior.time(behavior.trajbound==0),...
               tperfOutbound, behavior.time);

disp("Done with addToBehavior ... took " + num2str(toc) + " seconds");

% Let's add a call to compute local idphi
behavior.idphi = idphi(behavior.X, behavior.Y, 10); % 1/3 of a second integration window
