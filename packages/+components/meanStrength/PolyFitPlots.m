% all_times
load("PolyFit_MeanStrength.mat")

fig('mean strength all times, polyfit slopes, by pattern and dimension')

hpcSubset = PolyFit_MeanStrength.inx == "hpc-hpc";
hpcSubset = PolyFit_MeanStrength(hpcSubset,:);
pfcSubset = PolyFit_MeanStrength.inx == "hpc-pfc";
pfcSubset = PolyFit_MeanStrength(pfcSubset,:);
%%
% hpc-hpc
subplot(2,2,1)
 title("hpc-hpc, by pattern")
% compare the slope of all pats, high theta, delta, and ripple
mean_all_hpc = mean(str2double(hpcSubset.all),1);
mean_theta_hpc = mean(str2double(hpcSubset.theta),1);
mean_delta_hpc = mean(str2double(hpcSubset.delta),1);
mean_ripple_hpc = mean(str2double(hpcSubset.ripple),1);
std_all_hpc = std(str2double(hpcSubset.all));
std_theta_hpc = std(str2double(hpcSubset.theta));
std_delta_hpc = std(str2double(hpcSubset.delta));
std_ripple_hpc = std(str2double(hpcSubset.ripple));
errorbar([mean_all_hpc,mean_theta_hpc,mean_delta_hpc,mean_ripple_hpc],[std_all_hpc,std_theta_hpc,std_delta_hpc,std_ripple_hpc])
xticklabels(["all", "theta","delta","ripple"])

% hpc-hpc
subplot(2,2,2)
title("hpc-hpc, by leading dims")
% compare the slope of all pats, first 3 dimensions
mean_hpc_1 = mean(str2double(hpcSubset.dim1),1);
mean_hpc_2 = mean(str2double(hpcSubset.dim1),1);
mean_hpc_3 = mean(str2double(hpcSubset.dim1),1);

std_hpc_1 = std(str2double(hpcSubset.dim1));
std_hpc_2 = std(str2double(hpcSubset.dim2));
std_hpc_3 = std(str2double(hpcSubset.dim3));
bar([mean_all_hpc,mean_hpc_1,mean_hpc_2,mean_hpc_3]);
hold on
errorbar([mean_all_hpc,mean_hpc_1,mean_hpc_2,mean_hpc_3],[std_all_hpc,std_hpc_1,std_hpc_2,std_hpc_3])
xticklabels(["all", "1st dim","2nd dim","3rd dim"])

subplot(2,2,3)
title("hpc-pfc, by patterns")
% compare the slope of all pats, high theta, delta, and ripple
mean_all_pfc = mean(str2double(pfcSubset.all),1);
mean_theta = mean(str2double(pfcSubset.theta),1);
mean_delta_pfc = mean(str2double(pfcSubset.delta),1);
mean_ripple_pfc = mean(str2double(pfcSubset.ripple),1);
std_all_pfc = std(str2double(pfcSubset.all));
std_theta_pfc = std(str2double(pfcSubset.theta));
std_delta_pfc = std(str2double(pfcSubset.delta));
std_ripple_pfc = std(str2double(pfcSubset.ripple));
bar([mean_all_pfc,mean_theta,mean_delta_pfc,mean_ripple_pfc]);
hold on
errorbar([mean_all_pfc,mean_theta,mean_delta_pfc,mean_ripple_pfc],[std_all_pfc,std_theta_pfc,std_delta_pfc,std_ripple_pfc])
xticklabels(["all", "theta","delta","ripple"])


subplot(2,2,4)
 title("hpc-pfc, by leading dims")
% compare the slope of all pats, first 3 dimensions
mean_pfc_1 = mean(str2double(pfcSubset.dim1),1);
mean_pfc_2 = mean(str2double(pfcSubset.dim1),1);
mean_pfc_3 = mean(str2double(pfcSubset.dim1),1);

std_pfc_1 = std(str2double(pfcSubset.dim1));
std_pfc_2 = std(str2double(pfcSubset.dim2));
std_pfc_3 = std(str2double(pfcSubset.dim3));
errorbar([mean_all_pfc,mean_pfc_1,mean_pfc_2,mean_pfc_3],[std_all_pfc,std_pfc_1,std_pfc_2,std_pfc_3])
xticklabels(["all", "1st dim","2nd dim","3rd dim"])

