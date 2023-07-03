%% initialize data structures
%Patterns_AllAnimals = permute(Patterns,[2 1 3 4 5]);
%szPatterns = size(Patterns_AllAnimals);
%Patterns_AllAnimals = reshape(Patterns_AllAnimals, [szPatterns(1), prod(szPatterns(2:3)), szPatterns(4:5)]);

% numPairsPerPattern = (nSource+nTarget)^2 - nSource^2 - nTarget^2;
withpfc_pairs = cell(3,1);
withhpc_pairs = cell(3,1);

all_pairs_withpfc = [];
all_pairs_withhpc = [];

mean_corrwithpfc = [];
mean_corrwithhpc = [];

mean_withhpccorr_pattern = zeros(1,nPatterns);
mean_withpfccorr_pattern = zeros(1,nPatterns);

std_corrwithhpc = [];
std_corrwithpfc = [];


%% loop 
% and obtain the degree of co-firing within each and all patterns
nAnimalPartition = nAnimal * nPartition;

for i = progress(1:nPatterns)
    withhpc_pairs{i} = [];
    withpfc_pairs{i} = [];
    for p = progress(1:nAnimalPartition)
        % pick the source/target population
        curr_source = Patterns_AllAnimals(1,p,1,i).X_source';
        curr_HPC    = Patterns_AllAnimals(1,p,hpc,i).X_target';
        curr_PFC    = Patterns_AllAnimals(1,p,pfc,i).X_target';
        
        % calcualte co-firing of source and two target populations
        linearized_withhpc = calculatePatternCofiring(curr_source, curr_HPC);
        linearized_withpfc = calculatePatternCofiring(curr_source, curr_PFC);
        
        % append to all
        all_pairs_withhpc = [all_pairs_withhpc, linearized_withhpc];
        all_pairs_withpfc = [all_pairs_withpfc, linearized_withpfc];
        
        % append to pattern-specific
        withpfc_pairs{i} = [withpfc_pairs{i}, linearized_withpfc];
        withhpc_pairs{i} = [withhpc_pairs{i}, linearized_withhpc];
        
        mean_withhpccorr_pattern(i) = mean(linearized_withhpc(~isnan(linearized_withhpc)));
        mean_withpfccorr_pattern(i) = mean(linearized_withpfc(~isnan(linearized_withpfc)));
    end
    cur_mean_withhpc = mean(all_pairs_withhpc(~isnan(all_pairs_withhpc)));
    cur_std_withhpc = std(all_pairs_withhpc(~isnan(all_pairs_withhpc)));
    
    cur_mean_withpfc = mean(all_pairs_withpfc(~isnan(all_pairs_withpfc)));
    cur_std_withpfc = std(all_pairs_withpfc(~isnan(all_pairs_withpfc)));
    
    mean_corrwithhpc = [mean_corrwithhpc, cur_mean_withhpc];
    mean_corrwithpfc = [mean_corrwithpfc, cur_mean_withpfc];
    
    std_corrwithhpc = [std_corrwithhpc, cur_mean_withhpc];
    std_corrwithpfc = [std_corrwithpfc, cur_mean_withpfc];
end
