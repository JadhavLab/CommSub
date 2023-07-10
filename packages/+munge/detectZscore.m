function iszscore = detectZscore(data, p_struct)
% detectZscore - detect if data is zscored in a cheap easy way
% we either have rates/counts or zscored data ... so the test
% does not have to be sophisticated

assert(size(data,2) > size(data,1), 'Sample should be in columns');

iszscore = all(mean(data(1:3,:),2) < eps(class(data))*1e4);
if nargin == 1
    return
else
    if isfield(p_struct,'muFR')
        if ~iszscore
            warning('Data is not zscored, but muFR is provided.');
        end
        iszscore = true;
    else
        iszscore = false;
    end
end

end

