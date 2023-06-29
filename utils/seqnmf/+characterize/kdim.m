function dim = kdim(vartype)
%KDIM spits out the K dimension given a vartype, if it exists
switch vartype

case 'WH'
    dim = nan;
case 'W'
    dim = 2;
case 'H'
    dim = nan;
case ''
    dim = nan;

end
