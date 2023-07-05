function out = regionalexplained(Patterns, Option, varargin)
% this function calculates the regional firing prediction of all patterns
% across methods
% input: Patterns - a struct containing all the patterns
% output: out - a struct containing the regional firing prediction of all
% patterns across methods

ip = inputParser();
ip.addParameter('appendAttributes', {}, @(x) iscellstr(x) || isstring(x));
ip.parse(varargin{:});
Opt = ip.Results;

prophpc = struct();
proppfc = struct();
for prop = Opt.appendAttributes
    disp("init " + prop{1})
    prophpc.(prop{1}) = [];
    proppfc.(prop{1}) = [];
end

Patterns = squeeze(munge.reshapeByLabels(Patterns, 1, [Option.generateH],  'checksumSplitField', 'animal'));
Option   = munge.reshapeByLabels(Option, 1,   [Option.generateH], 'checksumSplitField', 'animal');
Patterns = permute(Patterns, [2 1 3 4 5]);
Patterns = reshape(Patterns, [size(Patterns,1), size(Patterns,2)*size(Patterns,3), size(Patterns,4), size(Patterns,5)]);

genHO = nd.fieldGet(Option, 'genH_name');
animalO = nd.fieldGet(Option, 'animal');
genHP = nd.fieldGet(Patterns, 'generateH');
animalP = nd.fieldGet(Patterns, 'animal');

hpc = Option(1).shortcut.HPC;
pfc = Option(1).shortcut.PFC;

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
patternVarExplained_hpc = zeros(nMethods, nPatterns, nAnimalPartition);
patternVarExplained_pfc = zeros(nMethods, nPatterns, nAnimalPartition);

genH = nd.fieldGet(Patterns, 'generateH');
animal = nd.fieldGet(Patterns, 'animal');
genH = squeeze(genH(:,1,1,:));
animal = squeeze(animal(:,:,1,:));
disp("genH: " + string(unique(genH)));

%% loop
for m = progress(1:nMethods, 'Title', char("method"))
    for i = 1:nPatterns
        for prop = Opt.appendAttributes
            prophpc(m,i).(prop{1}) = Patterns(m,1,hpc,i).(prop{1});
            proppfc(m,i).(prop{1}) = Patterns(m,1,pfc,i).(prop{1});
        end
        for p = progress(1:nAnimalPartition, 'Title', char("pattern " + i + " method " + m))

            % nCurrSource = Patterns(m,p,1,i).nSource;
            nCurrTarget = size(Patterns(m,p,2,i).X_target,1);
            % animalNo = floor(p/(nPartition+1)+1);
            % nCurrSource = nSource(animalNo);
            single_pred_with_hpc{m, i, p} = [];
            single_pred_with_pfc{m, i, p} = [];

            r_withhpc_patterns{m, i, p} = zeros(1,nCurrTarget);
            r_withpfc_patterns{m, i, p} = zeros(1,nCurrTarget);
            
            curr_source    = (Patterns(m,p,1,i).X_source)';
            curr_targethpc = (Patterns(m,p,hpc,i).X_target)';
            curr_targetpfc = (Patterns(m,p,pfc,i).X_target)';
            
            curr_B_hpc = Patterns(m,p,hpc,i).rankRegress.B;
            curr_B_pfc = Patterns(m,p,pfc,i).rankRegress.B;
            
            [patternhpc, meanhpc, ~] = ...
                plots.pred.var.explained(curr_source, curr_targethpc, curr_B_hpc);
            [patternpfc, meanpfc, ~] = ...
                plots.pred.var.explained(curr_source, curr_targetpfc, curr_B_pfc);
            
            patternVarExplained_hpc(m,i,p) = meanhpc;
            patternVarExplained_pfc(m,i,p) = meanpfc;
            
            r_withpfc_patterns{m, i, p} = patternpfc;
            r_withhpc_patterns{m, i, p} = patternhpc;

            % out.genH    = genH(i,m);
            % out.pattern = Option(m, i).patternNames(i);
        end
    end
end

out.r_withhpc_patterns      = r_withhpc_patterns; % regional firing prediction with hpc
out.r_withpfc_patterns      = r_withpfc_patterns; % regional firing prediction with pfc
out.single_pred_with_hpc    = single_pred_with_hpc; % single firing prediction with hpc
out.single_pred_with_pfc    = single_pred_with_pfc; % single firing prediction with pfc
out.patternVarExplained_hpc = patternVarExplained_hpc; % pattern variance explained with hpc
out.patternVarExplained_pfc = patternVarExplained_pfc; % variance explained with pfc
out.prophpc = prophpc;
out.proppfc = proppfc;

for prop = Opt.appendAttributes
    if length(unique([out.prophpc.(prop{1})])) == 1
        warning("prophpc." + prop{1} + " is constant")
    end
end

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
