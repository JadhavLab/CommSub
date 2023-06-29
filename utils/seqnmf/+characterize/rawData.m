function rawData(master, epoch_type)
% RAWDATA generates a plot of the raw data
%
if isequal(epoch_type,'run')
    yscale = size(master.data,1);
    imagesc(master.time, 1:yscale, master.data'); set(gca,'ydir','normal')
    cmocean('thermal')
    hold on
    plot(master.time, yscale*(master.trajdist + 1)/2, 'k', 'LineStyle', ':', 'LineWidth', 3)
    plot(master.time, smoothdata(yscale*master.vel./prctile(master.vel, 90)), 'LineStyle','--', 'Color', 'white')
else
    yscale = size(data.data,1);
    imagesc(master.time, 1:yscale, master.data'); set(gca,'ydir','normal')
    cmocean('thermal')
    hold on
end

xlim([6.3519    6.5982]*1000);
