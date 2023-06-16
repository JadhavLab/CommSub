load (animal+"pos01.mat");
pos = ndb.toNd(pos);
postime = [];
pos_x = [];
pos_y = [];

for i = 2:2:numel(pos)
    postime = [postime, pos(i).data(:,1)'];
    pos_x = [pos_x, pos(i).data(:,2)'];
    pos_y = [pos_y, pos(i).data(:,3)'];
end


running_spikeTimes = timeBinMidPoints(sessionTypePerBin == 1);
running_subspaces = total_subspaces(:,sessionTypePerBin == 1);
indicesToTranslate = BetweenTimes(postime, running_spikeTimes);
selectPosTime = postime(indicesToTranslate (~isnan(indicesToTranslate)));
selectPosX = pos_x(indicesToTranslate(~isnan(indicesToTranslate)));
selectPosY = pos_y(indicesToTranslate(~isnan(indicesToTranslate)));
%% calculate the location distribution of components
%
fig Location distrubtion; clf
tiledlayout(3,5); 
names = patternNames.replace('control','low');
for i = 1:3
    for j = 1:5
        nexttile; 
        log_track_distribution = signedlog(components.calculateLocationDistribution(selectPosX,selectPosY, 1, running_subspaces((i-1)*5 + j,:)));
        imagesc(log_track_distribution,'AlphaData',~isnan(log_track_distribution)); 
        crameri('vik', 'pivot', 0) % sets to vik colorscheme
        title(names(i) + " hpc-pfc {\it communication}" + newline +  "{\bf component_" + j + "}",'FontSize',8,'Interpreter','tex')
    end
end
sgtitle(animal+" top 5 components")

median_cenetered_runnning_subspaces = running_subspaces - median(running_subspaces,2);
fig Median-cenetedf location distrubtion; clf
tiledlayout(3,5); 
names = patternNames.replace('control','low');
for i = 1:3
    for j = 1:5
        nexttile; 
        log_track_distribution = signedlog(components.calculateLocationDistribution(selectPosX,selectPosY, 1, median_cenetered_runnning_subspaces((i-1)*5 + j,:)));
        imagesc(log_track_distribution,'AlphaData',~isnan(log_track_distribution)); 
        crameri('vik', 'pivot', 0) % sets to vik colorscheme
        title(names(i) + " hpc-pfc {\it communication}" + newline +  "{\bf component_" + j + "}",'FontSize',8,'Interpreter','tex')
    end
end
sgtitle(animal+" top 5 components")

TIME = 2;
zscore_runnning_subspaces = (running_subspaces - median(running_subspaces,TIME))./(std(running_subspaces,TIME));
fig Median-cenetedf location distrubtion; clf
tiledlayout(3,5); 
names = patternNames.replace('control','low');
for i = 1:3
    for j = 1:5
        nexttile; 
        log_track_distribution = signedlog(components.calculateLocationDistribution(selectPosX,selectPosY, 1, median_cenetered_runnning_subspaces((i-1)*5 + j,:)));
        imagesc(log_track_distribution,'AlphaData',~isnan(log_track_distribution)); 
        crameri('vik', 'pivot', 0) % sets to vik colorscheme
        title(names(i) + " hpc-pfc {\it communication}" + newline +  "{\bf component_" + j + "}",'FontSize',8,'Interpreter','tex')
    end
end
sgtitle(animal+" top 5 components")
