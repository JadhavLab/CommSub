function bar_comparison(means, B_high, B_low)

x = [1,2, 5,6];
colors = [cmocean('balance',4)];
temp = colors(1,:);
colors(1,:) = colors(2,:);
colors(2,:) = temp;
colors([1,4,5,8], :) =[];
figure;
for i = 1:4
    if i >= 3
        yyaxis right
    else
        yyaxis left
    end
    bar(x(i), means(i), 'facecolor', colors(i,:)); hold on;
    plot([x(i),x(i)], [B_low(i), B_high(i)], 'color', 0.5*colors(i,:),'linewidth', 2)
end

xticks(x)
xticklabels(["^{PFC}\delta\uparrow","^{PFC}\delta\downarrow","^{PFC}\delta\uparrow","^{PFC}\delta\downarrow"])
yyaxis right
ylabel('CA1 \theta')
yyaxis left
ylabel('WPLI \theta')
