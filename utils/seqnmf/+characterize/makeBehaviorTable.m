function [behWtab, varargout] = makeBehaviorTable(data)
% MAKEBEHAVIORTABLE creates a table of behavioral data for patterns K and times T

% Figure out my split datafields
dsources = ismember({'W_velocity','W_trajdist','W_aglindist'}, fields(data));

tabs     = {};
varnames = {};

%% Velocity
% --------
if dsources(1)
    [tabs{end+1}, varnames{end+1}] = makeSubtable(data, 'velocity');
end

%% Trajdist
% --------
if dsources(2)
    [tabs{end+1}, varnames{end+1}] = makeSubtable(data, 'trajdist');
end

%% Aglindist
% ---------
if dsources(3)
    [tabs{end+1}, varnames{end+1}] = makeSubtable(data, 'aglindist');
end

%% Return the separate table results
varargout = [{tabs}, {varnames}];
tabs = cellfun(@(x) unstackTable(x), tabs, 'UniformOutput', false);
try
    behWtab = table(tabs{:}, 'VariableNames', varnames);
catch ME
    behWtab = [];
end

% --------------------------------------------------------------------------
function [tab] = unstackTable(tab)

tab = removevars(tab,'P');
substacks      = {};
substack_names = {};
data = {'mu','stderr','d'};
for datum = data
    substacks{end+1} = unstack(removevars(tab, setdiff(data, datum)), datum,'PrettyL','NewDataVariableNames', sort(unique(tab.PrettyL)));
    %substacks{end}.measure = repmat(datum{1}, size(substacks{end}, 1), 1);
    substack_names{end+1} = datum{1};
end
tab = table(substacks{:}, 'VariableNames', data);

% --------------------------------------------------------------------------
function [tab, varname] = makeSubtable(data, field)
% Collapse onehot encoding

%tab{end+1} = mat2label(data.WH_velocity, data.orig_velocity, data.velocity);
wfield = ['W_' field];
[rawLabels, prettyLabels, rep] = labeler.label(data, field);
rawLabels = rawLabels(1:rep:end);
%dataSepByLabel = num2cell(data.W_velocity, [numel(rawLabels) 1 1]);
[P, K, T] = size(data.(wfield));
% Project repelem into new dimension
% and compute stats
% ----------------------------------
D = reshape(data.(wfield), [rep P/rep K T]);
mu     = squeeze(mean(D,1));
stderr = squeeze(std(D,1)/sqrt(rep));
P = P/rep;
assert(P == numel(rawLabels));
% Create table elements
% ----------------------------------
mu     = permute(mu,     [2 3 1]);
stderr = permute(stderr, [2 3 1]);
d = mu./stderr; 
[K, T, P] = ndgrid(1:K,1:T,1:P);
K = K(:);
T = T(:);
P = P(:);
R = string(rawLabels(P));
L = string(prettyLabels(P));
mu     = mu(:); % Ravel 3d (K X T X P) to 1d [PKT x 1]
stderr = stderr(:); % Ravel 3d (K X T X P) to 1d [PKT x 1]
d      = d(:); % Ravel 3d (K X T X P) to 1d [PKT x 1]
%dataSepByLabel = num2cell(dataSepByLabel', [size(dataSepByLabel,2) 1]);
tab = table(K, T, P, R, L, mu, stderr, d, 'VariableNames', ["K", "T", "P", "RawL", "PrettyL", "mu", "stderr", "d"]);
varname = field;
