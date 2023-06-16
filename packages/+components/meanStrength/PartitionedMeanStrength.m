%% poly-fits for changes in component strength mean, over all times

for j = 1:12
    for p = 1:Option.numPartition
        P_hpc_sv_all{j}{p} = polyfit(1:numel(hpc_m_all(p,j,:)),hpc_m_all(p,j,:),1);
        P_pfc_sv_all{j}{p} = polyfit(1:numel(pfc_m_all(p,j,:)),pfc_m_all(p,j,:),1);
        
        slope_hpc_sv_all(j,p) = P_hpc_sv_all{j}{p}(1);
        slope_pfc_sv_all(j,p) = P_pfc_sv_all{j}{p}(1);
    end
end
%%
%store into a table
animal = Option.animal;
animalName = repmat(animal, 2*Option.numPartition, 1);
inxName = [repmat("hpc-hpc", Option.numPartition, 1);repmat("hpc-pfc", Option.numPartition, 1)];
SlopeTable = array2table([animalName, inxName, [slope_hpc_sv_all' ; slope_pfc_sv_all']]);
SlopeTable.Properties.VariableNames(1:14) = ...
    {'animal','inx','all','theta','delta','ripple','low-theta','low-delta','low-ripple','dim1','dim2','dim3','dim4','dim5'};
load('PolyFit_MeanStrength.mat');
PolyFit_MeanStrength = [PolyFit_MeanStrength; SlopeTable];
save('PolyFit_MeanStrength','PolyFit_MeanStrength','-v7.3')
% Default heading for the columns will be A1, A2 and so on.
% You can assign the specific headings to your table in the following manner
%% poly-fits for changes in component strength mean, by behaviors