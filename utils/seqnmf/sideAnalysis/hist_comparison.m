function hist_comparison(B)


x = [1,2, 5,6];
colors = [cmocean('balance',4)];
temp = colors(1,:);
colors(1,:) = colors(2,:);
colors(2,:) = temp;

figure;
for i = 1:2:4
    if i == 1
        subplot(2,1,1)
    else
        subplot(2,1,2)
    end
    histogram(B{i}, 'facecolor', colors(i,:), 'facealpha', 0.25'); hold on
    histogram(B{i+1}, 'facecolor', colors(i+1,:), 'facealpha', 0.25')
end
