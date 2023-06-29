function [futurereward, pastreward, behavior] = futurePastReward(behavior, varargin)
% Adds future and past reward to behavior
%
% INPUT
% -----
% behavior : the table containing all information about the behavior
% group : the groups to compute the trial wise reward over 

ip=inputParser;
ip.addParameter('group',["animal","epoch"])
ip.addParameter('n',1)
ip.parse(varargin{:})
opt = ip.Results;

groupitems = arrayfun(@(x) behavior.(x), opt.group, 'UniformOutput', false);
G = findgroups(groupitems{:});

pastreward = zeros(numel(behavior.trial), opt.n);
for n = 1:opt.n
    for g = unique(G)'
        B = behavior(G==g,:);
        for trial = unique(B.trial)'
            if trial > n
                value = B.rew(B.trial == (trial-n));
                pastreward(behavior.trial == trial & G==g, n) = value(1) * ones(sum(B.trial==trial),1);
            else
                pastreward(behavior.trial == trial & G==g, n) = nan(sum(B.trial==trial),1);
            end
        end
    end
end
pastreward = pastreward(:,opt.n:-1:1);


futurereward = zeros(numel(behavior.trial), opt.n);
for n = 1:opt.n
for g = unique(G)'
    B = behavior(G==g,:);
    trials = unique(B.trial)';
    for trial = trials
        if trial < max(trials) - n + 1
            value = B.rew(B.trial == (trial+n));
            futurereward(behavior.trial == trial & G==g, n) = value(1) * ones(sum(B.trial==trial),1);
        else
            futurereward(behavior.trial == trial & G==g, n) = nan(sum(B.trial==trial),1);
        end
    end
end
end

behavior.futurereward = futurereward;
behavior.pastreward   = pastreward;
