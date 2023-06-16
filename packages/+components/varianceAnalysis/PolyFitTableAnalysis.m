% all_times
load("PolyFitTable_AllTimes.mat")

fig('all times, polyfit slopes, by pattern and dimension')

hpcSubset = PolyFitTable_AllTimes.inx == "hpc-hpc";
hpcSubset = PolyFitTable_AllTimes(hpcSubset,:);
pfcSubset = PolyFitTable_AllTimes.inx == "hpc-pfc";
pfcSubset = PolyFitTable_AllTimes(pfcSubset,:);
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





%%
% by_behaviors, such that the plots are sturctred the same as above but are
% for different behaviors
load("PolyFitTable_ByBeh.mat")




hpcSubset = PolyFitTable_ByBeh.inx == "hpc-hpc";
hpcSubset = PolyFitTable_ByBeh(hpcSubset,:);
pfcSubset = PolyFitTable_ByBeh.inx == "hpc-pfc";
pfcSubset = PolyFitTable_ByBeh(pfcSubset,:);
%%
behaviors = ["reward","error","inBoundChoice","outBoundChoice"];

for i = 1:4
    fig('polyfit slopes, by pattern and dimension when ' + behaviors(i)); clf;
    sgtitle(behaviors(i))
    beh_hpcSubset = hpcSubset.beh == behaviors(i);
    beh_pfcSubset = pfcSubset.beh == behaviors(i);
    beh_hpcSubset = hpcSubset(beh_hpcSubset,:);
        beh_pfcSubset = pfcSubset(beh_pfcSubset,:);


    
    % hpc-hpc
    subplot(2,2,1)
    
    % compare the slope of all pats, high theta, delta, and ripple
    mean_all_hpc = mean(str2double(beh_hpcSubset.all),1);
    mean_theta_hpc = mean(str2double(beh_hpcSubset.theta),1);
    mean_delta_hpc = mean(str2double(beh_hpcSubset.delta),1);
    mean_ripple_hpc = mean(str2double(beh_hpcSubset.ripple),1);
    std_all_hpc = std(str2double(beh_hpcSubset.all));
    std_theta_hpc = std(str2double(beh_hpcSubset.theta));
    std_delta_hpc = std(str2double(beh_hpcSubset.delta));
    std_ripple_hpc = std(str2double(beh_hpcSubset.ripple));
    bar([mean_all_hpc,mean_theta_hpc,mean_delta_hpc,mean_ripple_hpc]);
    hold on
    errorbar([mean_all_hpc,mean_theta_hpc,mean_delta_hpc,mean_ripple_hpc],[std_all_hpc,std_theta_hpc,std_delta_hpc,std_ripple_hpc])
    xticklabels(["all", "theta","delta","ripple"])
    title("hpc-hpc, by pattern")
    
    % hpc-hpc
    subplot(2,2,2)
    
    % compare the slope of all pats, first 3 dimensions
    mean_hpc_1 = mean(str2double(beh_hpcSubset.dim1),1);
    mean_hpc_2 = mean(str2double(beh_hpcSubset.dim1),1);
    mean_hpc_3 = mean(str2double(beh_hpcSubset.dim1),1);
    
    std_hpc_1 = std(str2double(beh_hpcSubset.dim1));
    std_hpc_2 = std(str2double(beh_hpcSubset.dim2));
    std_hpc_3 = std(str2double(beh_hpcSubset.dim3));
    bar([mean_all_hpc,mean_hpc_1,mean_hpc_2,mean_hpc_3]);
    hold on
    errorbar([mean_all_hpc,mean_hpc_1,mean_hpc_2,mean_hpc_3],[std_all_hpc,std_hpc_1,std_hpc_2,std_hpc_3])
    xticklabels(["all", "1st dim","2nd dim","3rd dim"])
    title("hpc-hpc, by leading dims")
    
    subplot(2,2,3)
    % compare the slope of all pats, high theta, delta, and ripple
    mean_all_pfc = mean(str2double(beh_pfcSubset.all),1);
    mean_theta = mean(str2double(beh_pfcSubset.theta),1);
    mean_delta_pfc = mean(str2double(beh_pfcSubset.delta),1);
    mean_ripple_pfc = mean(str2double(beh_pfcSubset.ripple),1);
    std_all_pfc = std(str2double(beh_pfcSubset.all));
    std_theta_pfc = std(str2double(beh_pfcSubset.theta));
    std_delta_pfc = std(str2double(beh_pfcSubset.delta));
    std_ripple_pfc = std(str2double(beh_pfcSubset.ripple));
    bar([mean_all_pfc,mean_theta,mean_delta_pfc,mean_ripple_pfc]);
    hold on
    errorbar([mean_all_pfc,mean_theta,mean_delta_pfc,mean_ripple_pfc],[std_all_pfc,std_theta_pfc,std_delta_pfc,std_ripple_pfc])
    xticklabels(["all", "theta","delta","ripple"])
    title("hpc-pfc, by patterns")
    
    
    subplot(2,2,4)
    
    % compare the slope of all pats, first 3 dimensions
    mean_pfc_1 = mean(str2double(beh_pfcSubset.dim1),1);
    mean_pfc_2 = mean(str2double(beh_pfcSubset.dim1),1);
    mean_pfc_3 = mean(str2double(beh_pfcSubset.dim1),1);
    
    std_pfc_1 = std(str2double(beh_pfcSubset.dim1));
    std_pfc_2 = std(str2double(beh_pfcSubset.dim2));
    std_pfc_3 = std(str2double(beh_pfcSubset.dim3));
    bar([mean_all_pfc,mean_pfc_1,mean_pfc_2,mean_pfc_3]);
    hold on
    errorbar([mean_all_pfc,mean_pfc_1,mean_pfc_2,mean_pfc_3],[std_all_pfc,std_pfc_1,std_pfc_2,std_pfc_3])
    xticklabels(["all", "1st dim","2nd dim","3rd dim"])
    title("hpc-pfc, by leading dims")
end
