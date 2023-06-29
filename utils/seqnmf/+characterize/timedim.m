function dim = timedim(vartype, field)
% TIMEDIM spits out the time dimension given a vartype

if isstring(vartype)
    vartype = char(vartype);
end
switch vartype
case ''
    if ~isempty(field)
        switch field
        case {'t','data'}
            dim = 2;
        otherwise
            dim = 1;
        end
    else
        dim = 2;
    end
case 'WH'
    dim = 2;
case 'W'
    dim = 3;
case 'H'
    dim = 2;
otherwise
    error(['Unrecognized vartype=' char(vartype)])
end
