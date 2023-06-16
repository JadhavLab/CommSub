function velocity_correlations = getPatternSpeedCorr(animal, H, Hvals, Htimes)


load(animal+"pos01.mat");
posdata = pos{1}{2}.data;
postime = unique(posdata(:,1));
posvel = posdata(:,end);
Hvals(:,3) = H(:,3);
Hvals(isnan(Hvals(:,3)),3) =0;
indices = BetweenTimes(Htimes, postime); % since H_times has the higher sampling rate, we want to translate the higher to the lower
indices = indices(~isnan(indices));

% posvel = interp1(postime, posvel, indices);
% Theta has positive speed correlation
velocity_correlations = corrcoef(...
    [posvel(1:numel(indices)), Hvals(indices,:)]);
velocity_correlations(isnan(velocity_correlations))=0;
end

