function [data, pre] = GetData(opt)
% Input
% -----
% variables in workspace at start of the script
%
% animal_list : cell of animals
debug=false;

% ----------------------------------------
% Iterate each animal and run seqnmf logic
% ----------------------------------------
P = []; % for aggreagting pos
C = []; % for aggreagting pos
T = []; % for lregion

cnt           = 0;
animal_cnt    = 0;
offset        = 0;
dayepoch_prev = [0 0];

for animal = opt.animal_list
    animal_cnt = animal_cnt + 1;


    % Determine the days and epochs for the animal
    animaldef_ = animaldef(animal{1});
    task = loaddatastruct(animaldef_{2:3}, 'task');
    search_results = cellfetch(task, 'environment');
    switch opt.epoch_type
    case 'run'
        epoch_filter = cellfun(@(x) isequal(x,'wtr1'), search_results.values);
    case 'sleep'
        epoch_filter = cellfun(@(x) isempty(x), search_results.values);
    case 'all'
        epoch_filter = ones(1,numel(search_results.values));
    otherwise
        error("Improper epoch type")
    end
    dayepochs = search_results.index(epoch_filter,:);
    if debug
        dayepochs
    end

    % ------------------------------
    % Load up all data for an animal
    % ------------------------------
    for dayepoch = dayepochs'

        day   = dayepoch(1);
        epoch = dayepoch(2);
        if sum(dayepoch - dayepoch_prev) < 0
            load(prev_avgeeg_file)
            avgeeg = ffend(avgeeg);
            if debug
                offset = offset  + avgeeg.endtime
            end
        end
        % ------------------------
        % -- load spectral info --
        % ------------------------
        if isequal(opt.epoch_type,'run')
            cgramfile = 'cgramcnew';
        else
            %cgramfile = 'cgramc';
            cgramfile = 'cgramcnew';
        end
        cgramc_file = [animal{1} cgramfile sprintf('-%02d-%02d.mat', dayepoch)];
        if ~exist(cgramc_file, 'file')
            disp('Continuing')
            continue
        else
            disp("Processing")
        end
        load(cgramc_file)
        cgramc = cgramcnew;
        data = ffend(cgramc);
        data.orig.t    = data.t;
        data.t         = data.t + offset;
        data.animalcnt = animal_cnt * ones(size(data.t));
        data.epoch     = epoch * ones(size(data.t));
        %data.t = data.t - avgeeg.starttime;
        % ------------------------
        % -- load position info --
        % ------------------------
        try
            C = [C; data];
            if ~isequal(opt.epoch_type, 'sleep')
                load([animal{1} 'pos' sprintf('%02d', day)])
                load([animal{1} 'lregion' sprintf('%02d', day)])
                p = pos{day}{epoch};
                p.data(:,1) = p.data(:,1) + offset;
                P = [P; repmat([animal_cnt, cnt],length(p.data),1), p.data];
                p = lregion{day}{epoch};
                p.time_orig = p.time; 
                p.time = p.time + offset;
                T = [T; repmat([animal_cnt,     cnt], length(p.time),1), p.time, p.TR, p.TB, p.lindist, p.lindist .* sign(p.TB(:,2)-0.5), p.tperf_timewise, p.time_orig]; % animal epoch time trial region rewarded trajbound lindist trajdist
            end
        catch ME
            warning(sprintf('Skipping animal %s, day %d epoch %d\n', animal{1}, dayepoch));
        end
        prev_avgeeg_file = [animal{1} 'avgeeg' sprintf('%02d-%02d.mat', dayepoch)];
        dayepoch_prev = dayepoch;
    end
end

% -------------------------------------
%     _                _               
%    / \   _ __   __ _| |_   _ _______ 
%   / _ \ | '_ \ / _` | | | | |_  / _ \
%  / ___ \| | | | (_| | | |_| |/ /  __/
% /_/   \_\_| |_|\__,_|_|\__, /___\___|
%                        |___/         
% -------------------------------------

if isequal(opt.epoch_type,'run')
    % -----------------------------------------
    % Convert pos and trial to common time
    % -----------------------------------------
    throwaway = ~ismember(P(:, 1:3), T(:, 1:3), 'rows');
    P(throwaway,:) = [];
    throwaway = ~ismember(T(:, 1:3), P(:, 1:3), 'rows');
    T(throwaway,:) = [];
    % Sort by animal epoch time
    P = sortrows(P, [1,2,3]);
    T = sortrows(T, [1,2,3]);
    % -------------------------------
    % Create convenience table arrays
    % -------------------------------
    assert(isequal(P(:,1:3),T(:,1:3)))
    % Behavior dataset
    P = num2cell(P,1);
    P = table(P{:}, ...
    'VariableNames', {'animal','epoch','time','x','y','dir','vel'});
    T = num2cell(T(:,4:end),1);
    T = table(T{:}, ...
    'VariableNames', {'trial','region','rew','traj', 'lindist', 'trajdist', 'tperf', 'orig_time'});
    behavior = [P,T];
    %data.T = [0 diff(data.t)];
    %data.T(data.T<0) = 0;
    %data.T = cumsum(data.T);
    %behavior.time(:) = [0; diff(behavior.time(:))];
    %behavior.time(behavior.time<0) = 0;
    %behavior.time = cumsum(behavior.time);
end

%--------------------------------------------------
%    / \   __ _  __ _ _ __ ___  __ _  __ _| |_ ___ 
%   / _ \ / _` |/ _` | '__/ _ \/ _` |/ _` | __/ _ \
%  / ___ \ (_| | (_| | | |  __/ (_| | (_| | ||  __/
% /_/   \_\__, |\__, |_|  \___|\__, |\__,_|\__\___|
%         |___/ |___/          |___/               
%--------------------------------------------------
% Time and frequency
data.t     = cat(2, C.t);
torig      = arrayfun(@(x) x.orig.t, C, 'UniformOutput', false);
%data.pre.t = cat(2, torig{:});
data.f     = C.f;
% Animal
% data.animal = cell(1,numel({C.t})); 
% cnt = 0;
% for t = {C.t}
%     cnt = cnt + 1;
%     data.animal{cnt} = cnt * ones(size(t{1}));
% end
data.animalcnt = cat(2, C.animalcnt);
data.animalnames = string(opt.animal_list);
data.epoch     = cat(2, C.epoch);
clear torig;

% Store potentially modified field list
data.opt = opt

Dt = diff(data.t); 
Dt(Dt<0)=0;
Dt(Dt>10)=0;
disp(['Total seconds of data: ' num2str(sum(Dt))])
clear Dt

% ADD HOTENCODED BEHAVIOR FIELDS
% ------------------------------
if any(contains(opt.fields, 'velocity'))
    behavior2data_inds = lookup(double(data.t), double(behavior.time));
    data.behavior2data_inds = behavior2data_inds;
    data.velocity = log(behavior.vel(behavior2data_inds));
    if opt.onehotencoding
        data.velocity_repelem = 3;
        data.orig_velocity    = data.velocity;
        data.velocity         = label2mat(data.velocity, 20, 'method', 'quantile')';
        data.velocity         = repelem(data.velocity, 1, data.velocity_repelem);
    end
end
if any(contains(opt.fields, 'trajdist'))
    behavior2data_inds = lookup(double(data.t), double(behavior.time));
    data.behavior2data_inds = behavior2data_inds;
    data.trajdist = behavior.trajdist(behavior2data_inds);
    if opt.onehotencoding
        data.trajdist_repelem = 3;
        data.orig_trajdist    = data.trajdist;
        data.trajdist         = label2mat(data.trajdist, 40)';
        data.trajdist         = repelem(data.trajdist, 1, data.trajdist_repelem);
    end
end
if any(contains(opt.fields, 'lindist'))
    behavior2data_inds = lookup(double(data.t), double(behavior.time));
    data.lindist = behavior.lindist(behavior2data_inds);
    if opt.onehotencoding
        data.lindist_repelem = 3;
        data.orig_lindist    = data.lindist;
        data.lindist         = label2mat(data.lindist, 20)';
        data.lindist         = repelem(data.lindist, 1, data.lindist_repelem);
    end
end
if any(contains(opt.fields, 'aglindist'))
    behavior2data_inds = lookup(double(data.t), double(behavior.time));
    data.behavior2data_inds = behavior2data_inds;
    data.aglindist = behavior.lindist(behavior2data_inds);
    inbound = behavior.traj(behavior2data_inds) == 1;
    data.aglindist(inbound) = -1*data.aglindist(inbound) + 1; % Map 1->0 traj to 0->1, like the outbound
    if opt.onehotencoding
        data.aglindist_repelem = 3;
        data.orig_aglindist    = data.aglindist;
        data.aglindist         = label2mat(data.aglindist, 10)';
        data.aglindist         = repelem(data.aglindist, 1, data.aglindist_repelem);
       
    end
end
if any(contains(opt.fields, 'trajbound'))
    behavior2data_inds = lookup(double(data.t), double(behavior.time));
    data.behavior2data_inds = behavior2data_inds;
    [~, ~, behavior] = characterize.futurePastTrajbound(behavior, 'n', 2);
    data.trajbound = [behavior.pasttrajbound(behavior2data_inds, :), behavior.rew(behavior2data_inds), behavior.futuretrajbound(behavior2data_inds, :)];
    if opt.onehotencoding
        data.trajbound_repelem = 1;
        data.orig_trajbound    = data.trajbound;
        %data.trajbound         = repelem(data.trajbound, 1, data.trajbound_repelem);
        data.trajbound         = arrayfun(@(col) label2mat(data.trajbound(:,col)+1',2,'method','raw'), ...
                                1:size(data.trajbound,2), 'UniformOutput', false);
        data.trajbound         = cat(1, data.trajbound{:})';
    end
end
if any(contains(opt.fields, 'reward'))
    behavior2data_inds = lookup(double(data.t), double(behavior.time));
    data.behavior2data_inds = behavior2data_inds;
    [~, ~, behavior] = characterize.futurePastReward(behavior, 'n', 2);
    data.reward = [behavior.pastreward(behavior2data_inds,:), behavior.rew(behavior2data_inds), behavior.futurereward(behavior2data_inds,:)];
    if opt.onehotencoding
        data.reward_repelem = 1;
        data.orig_reward    = data.reward;
        %data.reward         = repelem(data.reward, 1, data.reward_repelem);
        data.reward         = arrayfun(@(col) label2mat(data.reward(:,col)+1',2,'method','raw'), ...
                                1:size(data.reward,2), 'UniformOutput', false);
        data.reward         = cat(1, data.reward{:})';
    end
end
% Other opt.fields
for field = opt.fields
    if isfield(C,field{1})
        data.(field{1}) = cat(1, C.(field{1}));
    end
end
F = opt.fields; % tmp opt.fields for later, because phi maps to phi_cos and phi_sin
if any(contains(opt.fields, 'phi'))
    data.phi_sin = sin(data.phi);
    data.phi_cos = cos(data.phi);
    F{cellfun(@(x) isequal(x, 'phi'), F)} = 'phi_sin';
    F = [F 'phi_cos'];
end

% Pack opt.fields
% -----------
fieldpack_tmp = [opt.fieldpack_kws, 'groups', data.animalcnt];
data.data = seq.packfields(data, F, fieldpack_tmp{:})';

% Store copy of data before subsampling
% -------------------------------------
if opt.storeDatBeforeSubsamp && opt.subsample > 0
    for field = union(data.opt.fields, {'S1','S2','plv','wpli','C','trajdist','f','t','data'})
        if ~ismember(field{1}, behavior.Properties.VariableNames)
            pre.(field{1}) = data.(field{1});
        end
    end
else
    pre = [];
end

% SUBSAMPLE
% ---------
disp('Before')
size(data.data)
data = seq.subsampleByGroup(data, opt.fields, opt.subsample, 'storeOrig', opt.storeDatBeforeSubsamp);
disp('After')
size(data.data)
assert(numel(data) > 0);


clear P T cgramcnew cgramc behavior2data_inds pos task
