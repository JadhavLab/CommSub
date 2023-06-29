function gstruct = findgroups(tab, conditionLabels)
% Accepts either a table or a struct that looks like a table

if isstruct(tab)
    tab = struct2table(tab);
end

conditionLabels = string(conditionLabels);

conditionals              = cell(1,numel(conditionLabels));
conditions                = num2cell(table2array(tab(:, conditionLabels)),1);
addressConditionals       = num2cell(table2array(tab(:, conditionLabels)),1);
[groups, conditionals{:}] = findgroups(conditions{:});
uGroups                   = unique(groups);

gstruct.conditionLabels = conditionLabels;
gstruct.uGroups         = uGroups;
gstruct.time.groups     = groups;
gstruct.group.values    = conditionals;

addressConditionals = {};
for l = 1:numel(conditionLabels)
    label = conditionLabels(l);
    gstruct.time.field.(label) = conditionals{l}(groups);
    if isnumeric(conditionals{l})
        isAPositiveInt = isequal(round(conditionals{l}), conditionals{l});
        if isAPositiveInt
            isAPositiveInt = isAPositiveInt && all(conditionals{l} > 0, 'all');
        end
    else
        isAPositiveInt = false;
    end
    if isAPositiveInt
        addressConditionals{l} = gather(conditionals{l});
    else
        addressConditionals{l} = int32(categorical(conditionals{l}));
    end
    gstruct.group.field.(label) = nan(1, numel(uGroups));
    for g = 1:numel(uGroups)
        gstruct.group.field.(label)(g) = conditionals{l}(g);
    end
end
gstruct.group.address           = addressConditionals;
gstruct.group.addressByGroupNum = cat(2, addressConditionals{:});
gstruct.group.addressByGroupNum = num2cell(gstruct.group.addressByGroupNum, 2);

for g = uGroups'
    tmp = ["$"] + conditionLabels + " == " + string(gstruct.group.addressByGroupNum{g});
    gstruct.evalstring(g) = join(tmp, " & ");
end
