function out = regionalexplained(Patterns, Option)
% this function calculates the regional firing prediction of all patterns
% across methods
% input: Patterns - a struct containing all the patterns
% output: out - a struct containing the regional firing prediction of all
% patterns across methods

Patterns = munge.reshapeByLabels(Patterns, 1, [Option.generateH]);
szPatterns       = size(Patterns);
nMethods         = szPatterns(1);
nAnimalPartition = szPatterns(2);
nPatterns        = szPatterns(end);

% this script calculates the regional firing predicition of all patterns
% across methods
%% initialize data sturctures
r_withhpc_patterns   = cell(nMethods, nPatterns, nAnimalPartition);
r_withpfc_patterns   = cell(nMethods, nPatterns, nAnimalPartition);
single_pred_with_hpc = cell(nMethods, nPatterns, nAnimalPartition);
single_pred_with_pfc = cell(nMethods, nPatterns, nAnimalPartition);

%% loop
for m = 1:nMethods
    for i = 1:nPatterns
        for p = 1:nAnimalPartition
            % nCurrSource = Patterns(m,p,1,i).nSource;
            nCurrTarget = Patterns(m,p,1,i).nTarget;
            % animalNo = floor(p/(nPartition+1)+1);
            % nCurrSource = nSource(animalNo);
            single_pred_with_hpc{m}{i}{p} = [];
            single_pred_with_pfc{m}{i}{p} = [];

            r_withhpc_patterns{m}{i}{p} = zeros(1,nCurrTarget);
            r_withpfc_patterns{m}{i}{p} = zeros(1,nCurrTarget);
            
            curr_source = (Patterns(m,p,1,i).X_source)';
            curr_targethpc = (Patterns(m,p,hpc,i).X_target)';
            curr_targetpfc = (Patterns(m,p,pfc,i).X_target)';
            
            curr_B_hpc = Patterns(m,p,hpc,i).rankRegress.B;
            curr_B_pfc = Patterns(m,p,pfc,i).rankRegress.B;
            
            [patternhpc, meanhpc, ~] = calculateVarianceExplained...
                                      (curr_source, curr_targethpc, curr_B_hpc);
            [patternpfc, meanpfc, ~] = calculateVarianceExplained...
                                      (curr_source, curr_targetpfc, curr_B_pfc);
            
            patternVarExplained_hpc(m,i,p) = meanhpc;
            patternVarExplained_pfc(m,i,p) = meanpfc;
            
            r_withpfc_patterns{m}{i}{p} = patternpfc;
         
            r_withhpc_patterns{m}{i}{p} = patternhpc;
            

        end
    end
end

out.r_withhpc_patterns = r_withhpc_patterns;
out.r_withpfc_patterns = r_withpfc_patterns;
out.single_pred_with_hpc = single_pred_with_hpc;
out.single_pred_with_pfc = single_pred_with_pfc;
out.patternVarExplained_hpc = patternVarExplained_hpc;
out.patternVarExplained_pfc = patternVarExplained_pfc;

end

%% graveyard -- used to be inside loop
%             for k = 1:nCurrSource
%                 curr_singleB = Patterns...
%                                (m,p,1,i).rankRegress.singlesource_B{k};
%                 if ~isempty(curr_singleB)
%                      curr_singlesource = curr_source(:,k);
%                     [single_hpc,~] = calculatePredictionPerformance...
%                                               (curr_singlesource, curr_targethpc, curr_singleB);
%                     [single_pfc,~] = calculatePredictionPerformance...
%                                               (curr_singlesource, curr_targetpfc, curr_singleB);
%                     single_pred_with_hpc{m}{i}{p} = ...
%                                      [single_pred_with_hpc{m}{i}{p}, single_hpc];
%                     single_pred_with_pfc{m}{i}{p} =...
%                                      [single_pred_with_pfc{m}{i}{p}, single_pfc];
%                 end
%             end
