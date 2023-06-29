function [futuretrajbound, pasttrajbound, behavior] = futurePastTrajbound(behavior, varargin)
% Adds future and past trajbound to behavior
%
% INPUT
% -----
% behavior : the table containing all information about the behavior
% group : the groups to compute the trial wise trajbound over 

ip=inputParser;
ip.addParameter('group',["animal","epoch"])
ip.addParameter('n',1)
ip.parse(varargin{:})
opt = ip.Results;

groupitems = arrayfun(@(x) behavior.(x), opt.group, 'UniformOutput', false);
G = findgroups(groupitems{:});

pasttrajbound = zeros(numel(behavior.trial), opt.n);
for n = 1:opt.n
    for g = unique(G)'
        B = behavior(G==g,:);
        for trial = unique(B.trial)'
            if trial > n
                value = B.traj(B.trial == (trial-n));
                pasttrajbound(behavior.trial == trial & G==g, n) = value(1) * ones(sum(B.trial==trial),1);
            else
                pasttrajbound(behavior.trial == trial & G==g, n) = nan(sum(B.trial==trial),1);
            end
        end
    end
end
pasttrajbound = pasttrajbound(:,opt.n:-1:1);


futuretrajbound = zeros(numel(behavior.trial), opt.n);
for n = 1:opt.n
for g = unique(G)'
    B = behavior(G==g,:);
    trials = unique(B.trial)';
    for trial = trials
        if trial < max(trials) - n + 1
            value = B.traj(B.trial == (trial+n));
            futuretrajbound(behavior.trial == trial & G==g, n) = value(1) * ones(sum(B.trial==trial),1);
        else
            futuretrajbound(behavior.trial == trial & G==g, n) = nan(sum(B.trial==trial),1);
        end
    end
end
end

behavior.futuretrajbound = futuretrajbound;
behavior.pasttrajbound   = pasttrajbound;
