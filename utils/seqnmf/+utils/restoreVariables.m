function restoreVariables(S, where)
if nargin == 1
where = 'base';
end
disp('restoring variables')
for field = fieldnames(S)'
if ~isa(S.(field{1}),'matlab.ui.figure')
    try
        assignin(where,field{1},S.(field{1}))
    catch MatlabException
        warning('cannot load %s', field{1})
    end
end
end
close all % this is a hack. since I save figures in the matfile, they will appear whenever restoring.
disp('restored.')