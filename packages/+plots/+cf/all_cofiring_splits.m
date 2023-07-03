M=matfile('~/savestate.mat')

genH = shortcut.generateH(string(Option.generateH));
patternSymbols = shortcut.patternSymbols(patternNames,2);

spikeRateMatrix = M.spikeRateMatrix;
areasPerNeuron  = M.areaPerNeuron;

H = areasPerNeuron == "CA1";
P = areasPerNeuron == "PFC";
HH = {find(H),find(H)};
PP = {find(P),find(P)};
HP = {find(H),find(P)};

ravel = @(x) x(:);
R = corrcoef(spikeRateMatrix');
for i = 1:size(R,1); R(i,i) = nan; end
same = [ravel(R(HH{:})); ravel(R(PP{:}))];
same_hh = [ravel(R(HH{:}))];
same_pp = [ravel(R(PP{:}))];
different = [ravel(R(HP{:}))];

figc checksum
[N,edges] = histcounts(same);
N = N./sum(N);
bar(edges(1:end-1),N);
hold on
[N,edges] = histcounts(different);
N = N./sum(N);
bar(edges(1:end-1),N);
legend('Same','Different')

nexttile
n=50;
[N,edges] = histcounts(same_hh,n);
N = N./sum(N);
bar(edges(1:end-1),N);
hold on
[N,edges] = histcounts(same_pp,edges);
N = N./sum(N);
bar(edges(1:end-1),N);
hold on
[N,edges] = histcounts(different,edges);
N = N./sum(N);
bar(edges(1:end-1),N);
legend('Same_hh','Same_pp','Different')
s(1) = skewness(same_hh);
s(2) = skewness(same_pp);
s(3) = skewness(different);
q(1) = quantile(same_hh,0.60);
q(2) = quantile(same_pp,0.60);
q(3) = quantile(different,0.60);
message = ["skewness_{","quant^{60}_{"];
message = repmat(message,3,1) + repmat(["hh","pp","hp"]',1,2) + repmat("} = ",3,2) + [s',q'];
title(message)
set(gca,'yscale','log');
X_hpc = M.X_hpc;
X_pfc = M.X_pfc;

figc("Also true within rhythms?");
tiledlayout(2,6)
for i = progress(1:2)
    for j = progress(1:6)

        H = 1:size(X_hpc{i}{j},1);
        P = H(end) + 1:size(X_pfc{i}{j},1);
        HH = {H,H};
        PP = {P,P};
        HP = {H,P};        

        R = corrcoef([X_hpc{i}{j}', X_pfc{i}{j}']);
        for k = 1:size(R,1); R(k,k) = nan; end
        same      = [ravel(R(HH{:})); ravel(R(PP{:}))];
        same_hh   = [ravel(R(HH{:}))];
        same_pp   = [ravel(R(PP{:}))];
        different = [ravel(R(HP{:}))];

        nexttile
        [N,edges] = histcounts(same);
        edges = linspace(-0.1, 1,  round(1/median(diff(edges))));
        [N,edges] = histcounts(same, edges);
        N = N./sum(N);
        b= bar(edges(1:end-1),N);
        set(b,'FaceAlpha',0.33);
        hold on
        [N,edges] = histcounts(different,edges);
        N = N./sum(N);
        b=bar(edges(1:end-1),N);
        set(b,'FaceAlpha',0.33);
        legend('HPC-HPC/PFC-PFC','HPC-PFC')

        title(["GenH = " + genH(i), "Pat = " + patternSymbols(j)])

    end
end


figc("Also true within rhythms? split same");
tiledlayout(2,6)
for i = progress(1:2)
    for j = progress(1:6)

        H = 1:size(X_hpc{i}{j},1);
        P = H(end) + 1:size(X_pfc{i}{j},1);
        HH = {H, H};
        PP = {P, P};
        HP = {H, P};        

        Y = [X_hpc{i}{j}', X_pfc{i}{j}'];
        R = corrcoef(Y);
        for k = 1:size(R,1); R(k,k) = nan; end
        same      = [ravel(R(HH{:})); ravel(R(PP{:}))];
        same_hh   = [ravel(R(HH{:}))];
        same_pp   = [ravel(R(PP{:}))];
        different = [ravel(R(HP{:}))];

        nexttile
        [N,edges] = histcounts(same);
        edges = linspace(-0.1, 1,  round(1/median(diff(edges))));
        [N,edges] = histcounts(same_hh, edges);
        N = N./sum(N);
        b= bar(edges(1:end-1),N);
        set(b,'FaceAlpha',0.33);
        hold on
        [N,edges] = histcounts(same_pp,edges);
        N = N./sum(N);
        b=bar(edges(1:end-1),N);
        set(b,'FaceAlpha',0.33);
        hold on
        [N,edges] = histcounts(different,edges);
        N = N./sum(N);
        b=bar(edges(1:end-1),N);
        set(b,'FaceAlpha',0.33);
        legend('HPC-HPC','PFC-PFC','HPC-PFC')
        s(1) = skewness(same_hh);
        s(2) = skewness(same_pp);
        s(3) = skewness(different);
        q(1) = quantile(same_hh,0.85);
        q(2) = quantile(same_pp,0.85);
        q(3) = quantile(different,0.85);
        message = ["skewness_{","quant^{85}_{"];
        message = repmat(message,3,1) + repmat(["hh","pp","hp"]',1,2) + repmat("} = ",3,2) + [s',q'];
        title([message;"GenH = " + genH(i), "Pat = " + patternSymbols(j)])

    end
end
