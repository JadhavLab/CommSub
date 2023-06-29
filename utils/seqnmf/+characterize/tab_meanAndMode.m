function varargout = tab_meanAndMode(varargin)
% Carries out mean and mode automatically on a table object
%

varname = {};
tabs    = {};

T = varargin;

for col = 1:size(T, 2)
    
    if istable(T)
        avar = table2array(T(:,col));
    elseif iscell(T)
        avar = T{col};
    else
        avar = T(:,col);
    end

    if isnumeric(avar)
        tabs{end+1} = mean(avar, 1);
    else
        tabs{end+1} = mode(categorical(avar), 1);
    end

    if istable(T)
        varname{end+1} = T.Properties.VariableNames{col};
    end

end

keyboard
if istable(T)
    varargout{1} = table(tabs{:}, 'VariableNames', varname);
else
    varargout = tabs;
end

