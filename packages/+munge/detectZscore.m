function iszscore = detectZscore(data)
% detectZscore - detect if data is zscored in a cheap easy way
% we either have rates/counts or zscored data ... so the test
% does not have to be sophisticated

assert(size(data,2) > size(data,1), 'Sample should be in columns');

iszscore = all(mean(data(1:3,:),1) < eps(class(data))*10);

end

