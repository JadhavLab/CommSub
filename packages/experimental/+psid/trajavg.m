function avg= plotavg(b, resolution, varargin)
% Takes an behavior table annontated with X_pred

% Optionals
ip = inputParser;
ip.addParameter('propnames', ["epoch", "rewarded", "trajbound", "leftright", "ldistbin"]);
ip.addParameter('func', @nanmean);
ip.addParameter('propadd', ["field","col"]);
ip.addParameter('ploton', true);
ip.parse(varargin{:})
Opt = ip.Results;

% Determine props to avg
propadd   = string(Opt.propadd(ismember(Opt.propadd, fieldnames(b))));
propnames = string(Opt.propnames(ismember(Opt.propnames, fieldnames(b))));
propnames = [propnames, propadd];
propnames = [propnames, "ldistbin"];

% Discretize axis
discaxis = linspace(0, 1, resolution+1);
b.ldistbin = discretize(b.lindist, discaxis);

% Properties to split by
props = cell(numel(propnames),1);
propcell  = arrayfun(@(field) b.(field), propnames, 'UniformOutput', false);
[groups, props{:}] = findgroups(propcell{:});
uGroups = unique(groups);

% Find average trajects
avg = [];
for group = progress(uGroups', 'Title','Property averaging')

    filt = group == groups;
    pcnt = 0;
    for prop = propnames
        pcnt = pcnt + 1;
        avg(group).(prop) = props{pcnt}(group);
    end

    avg(group).val = Opt.func(b.val(filt));
    avg(group).N = sum(filt);

end

avg = struct2table(avg);

if Opt.ploton
    psid.plottrajavg(avg);
end
