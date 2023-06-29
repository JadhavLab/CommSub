function [windowsC, windowsS1, meanC, meanS1, center, windowsCtheta, windowsS1theta, meanCtheta, meanS1theta] = average_peaks(C)

delta = C.f >= 0.5 & C.f < 5;

before=10/C.params.movingwin(2);
after=10/C.params.movingwin(2);

mean_delta = mean(C.S2(:,delta),2);
delta_peak = 40;
[dpeaks, dlocs] = findpeaks(80*log10(mean_delta+10)/max(log10(mean_delta+10)), 'minpeakheight', delta_peak);

filt = dlocs > before & (dlocs <  size(C.S1, 1) - after - 1);
dpeaks = dpeaks(filt);
dlocs = dlocs(filt);


windowsC = cell(1,numel(dpeaks));
windowsS1 = cell(1,numel(dpeaks));
for i = 1:numel(dpeaks)

    peak = dpeaks(i);
    loc = dlocs(i);
    
    windowsC{i}  = C.wpli(loc-before:loc+after,:);
    windowsS1{i} = C.S1(loc-before:loc+after,:);

end

windowsC = cat(3, windowsC{:});
windowsS1 = cat(3, windowsS1{:});
meanC = mean(windowsC, 3);
meanS1 = mean(windowsS1, 3);

theta = C.f >= 5.5 & C.f < 13;

windowsCtheta = squeeze(mean(windowsC(:,theta,:),2));
windowsS1theta = squeeze(mean(windowsS1(:,theta,:),2));
meanCtheta = mean(windowsCtheta, 2);
meanS1theta = mean(windowsS1theta, 2);



center = false(before+after+1,1)';
center(before+1) = true;
