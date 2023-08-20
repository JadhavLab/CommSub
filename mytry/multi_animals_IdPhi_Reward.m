warning("off","all")

animal_lst = ["JS13", "JS15", "JS17", "JS21", "ZT2"];

for i = 4:5
disp("!!!!!!!!!!!!!!!!")
disp("The "+i+" Animal")
disp("!!!!!!!!!!!!!!!!")

animal = animal_lst(i);
[Option, Spk, Patterns, Patterns_overall] = analysis.setup_idphi_reward(animal);

a1 = size(Spk.hpc.T{1,1});
a2 = size(Spk.hpc.T{1,2});
a3 = size(Spk.hpc.T{1,3});
a4 = size(Spk.hpc.T{1,4});
trial = [a1(3),a2(3),a3(3),a4(3)];

disp("------------------------")
%% 
disp("        jpecc           ")
disp("------------------------")

Option.analysis.JPECC = 1;


J_Patterns.(animal) = analysis.JPECC(Patterns, Option, 1, trial);
J_Patterns2.(animal) = analysis.JPECC(Patterns, Option, 2, trial);
J_Patterns3.(animal) = analysis.JPECC(Patterns, Option, 3, trial);
J_Patterns4.(animal) = analysis.JPECC(Patterns, Option, 4, trial);

end




%%%%%%Gaussian Filter%%%%%%
sigma = 0.7;
filter_size = 2;
[x, y] = meshgrid(-filter_size:filter_size, -filter_size:filter_size);
h = exp(-(x.^2 + y.^2) / (2 * sigma^2));
h = h / sum(h(:));  % normalize the filter

%%%%%%Set for Plotting%%%%%%%%
stackedMatrices = cell(4, 40);
meanMatrices = cell(1,4);
for j = 1:4
    P = J_Patterns.(animal_lst(j));
for i = 1:40
    x = nd.fieldGet(P(i).jpecc, "val1");
    stackedMatrices{j,i} = x;
end
stackedMatrice = cat(40, stackedMatrices{j,:});
meanMatrices{j} = nanmean(stackedMatrice, 40);
end

stackedMatrices2 = cell(4, 40);
meanMatrices2 = cell(1,4);
for j = 1:4
    P = J_Patterns2.(animal_lst(j));
for i = 1:40
    x = nd.fieldGet(P(i).jpecc, "val1");
    stackedMatrices2{j,i} = x;
end
stackedMatrice2 = cat(40, stackedMatrices2{j,:});
meanMatrices2{j} = nanmean(stackedMatrice2, 40);
end

stackedMatrices3 = cell(4, 40);
meanMatrices3 = cell(1,4);
for j = 1:4
    P = J_Patterns3.(animal_lst(j));
for i = 1:40
    x = nd.fieldGet(P(i).jpecc, "val1");
    stackedMatrices3{j,i} = x;
end
stackedMatrice3 = cat(40, stackedMatrices3{j,:});
meanMatrices3{j} = nanmean(stackedMatrice3, 40);
end

stackedMatrices4 = cell(4, 40);
meanMatrices4 = cell(1,4);
for j = 1:4
    P = J_Patterns4.(animal_lst(j));
for i = 1:40
    x = nd.fieldGet(P(i).jpecc, "val1");
    stackedMatrices4{j,i} = x;
end
stackedMatrice4 = cat(40, stackedMatrices4{j,:});
meanMatrices4{j} = nanmean(stackedMatrice4, 40);
end

%%%%%%%%%%Plot IdPhi_High%%%%%%%%
for i = 1:4
Matrix = meanMatrices3{i};
%Matrix = conv2(Matrix, h, 'same');

% heatmap
figure;
t = pcolor(Matrix);
%shading interp;

caxis([0 0.5]);
cmocean('haline')
colorbar

%{
%caxis([-0.4 0.4]);
colormap(cmap);
colorbar;
%}

set(t, 'EdgeColor', 'none');
line([1 60], [1 60], 'Color', [1,1,1], 'LineWidth', 0.8,"LineStyle", "--");

xlabel('PFC');
ylabel('HPC');
title('IdPhi-High-' + animal_lst(i));
saveas(gcf, "F:\ComSub\Figures\multi_animals\idphi and reward\IdPhi-High-" + animal_lst(i) + ".png");
end

%%%%%%%%%%Plot IdPhi_Low%%%%%%%%
for i = 1:4
Matrix = meanMatrices4{i};
%Matrix = conv2(Matrix, h, 'same');

% heatmap
figure;
t = pcolor(Matrix);
%shading interp;

caxis([0 0.5]);
cmocean('haline')
colorbar

%{
%caxis([-0.4 0.4]);
colormap(cmap);
colorbar;
%}

set(t, 'EdgeColor', 'none');
line([1 60], [1 60], 'Color', [1,1,1], 'LineWidth', 0.8,"LineStyle", "--");

xlabel('PFC');
ylabel('HPC');
title('IdPhi-Low-' + animal_lst(i));
saveas(gcf, "F:\ComSub\Figures\multi_animals\idphi and reward\IdPhi-Low-" + animal_lst(i) + ".png");
end

%%%%%%%%%%Plot Reward_On%%%%%%%%
for i = 1:4
Matrix = meanMatrices{i};
%Matrix = conv2(Matrix, h, 'same');

% heatmap
figure;
t = pcolor(Matrix);
%shading interp;

caxis([-0.1 0.3]);
cmocean('haline')
colorbar

%{
%caxis([-0.4 0.4]);
colormap(cmap);
colorbar;
%}

set(t, 'EdgeColor', 'none');
line([1 60], [1 60], 'Color', [1,1,1], 'LineWidth', 0.8,"LineStyle", "--");

xlabel('PFC');
ylabel('HPC');
title('Reward-On-' + animal_lst(i));
saveas(gcf, "F:\ComSub\Figures\multi_animals\idphi and reward\Reward-On-" + animal_lst(i) + ".png");
end

%%%%%%%%%%Plot Reward_Off%%%%%%%%
for i = 1:4
Matrix = meanMatrices2{i};
%Matrix = conv2(Matrix, h, 'same');

% heatmap
figure;
t = pcolor(Matrix);
%shading interp;

caxis([-0.1 0.3]);
cmocean('haline')
colorbar

%{
%caxis([-0.4 0.4]);
colormap(cmap);
colorbar;
%}

set(t, 'EdgeColor', 'none');
line([1 60], [1 60], 'Color', [1,1,1], 'LineWidth', 0.8,"LineStyle", "--");

xlabel('PFC');
ylabel('HPC');
title('Reward-Off-' + animal_lst(i));
saveas(gcf, "F:\ComSub\Figures\multi_animals\idphi and reward\Reward-Off-" + animal_lst(i) + ".png");
end