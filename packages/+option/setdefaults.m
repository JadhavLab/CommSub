function Option = setdefaults(Option, name)


if nargin < 2
    name = 'TheScript';
end

Default = option.defaults(name)
for field = string(fieldnames(Default))'
    if ~isfield(Option, field)
        Option.(field) = Default.(field);
    end
end

if contains(Option.generateH, "fromRipTimes")
    Option.generateFromRipTimes = true;
else
    Option.generateFromRipTimes = false;
end

Option = option.postamble(Option);
