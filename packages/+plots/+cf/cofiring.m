function out = calculatePatternCofiring(Patterns, Option, varargin)

ip = inputParser();
ip.addParameter('appendAttributes', {}, @(x) iscellstr(x) || isstring(x));
ip.parse(varargin{:});
Opt = ip.Results;

%% initialize data structures
%Patterns = permute(Patterns,[2 1 3 4 5]);
%szPatterns = size(Patterns);
%Patterns = reshape(Patterns, [szPatterns(1), prod(szPatterns(2:3)), szPatterns(4:5)]);
% numPairsPerPattern = (nSource+nTarget)^2 - nSource^2 - nTarget^2;

Patterns = munge.resh12(Patterns);

short = Option(1).shortcut;
hpc = short.HPC;
pfc = short.PFC;

nPatterns = size(Patterns,3);
nAnimalPartition = size(Patterns,1);

%% allocate
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

prophpc = struct();
proppfc = struct();
for prop = Opt.appendAttributes
    disp("init " + prop{1})
    prophpc.(prop{1}) = [];
    proppfc.(prop{1}) = [];
end


%% loop 
% and obtain the degree of co-firing within each and all patterns
for i = progress(1:nPatterns)
    withhpc_pairs{i} = [];
    withpfc_pairs{i} = [];
    for p = progress(1:nAnimalPartition)
        % pick the source/target population
        curr_source = Patterns(p,1,i).X_source';
        curr_HPC    = Patterns(p,hpc,i).X_target';
        curr_PFC    = Patterns(p,pfc,i).X_target';

        for prop = Opt.appendAttributes
            prophpc.(prop{1}) = [prophpc.(prop{1}), Patterns(p,hpc,i).(prop{1})];
            proppfc.(prop{1}) = [proppfc.(prop{1}), Patterns(p,pfc,i).(prop{1})];
        end
        
        % calcualte co-firing of source and two target populations
        linearized_withhpc = plots.cf.calculatePatternCofiring(curr_source, curr_HPC);
        linearized_withpfc = plots.cf.calculatePatternCofiring(curr_source, curr_PFC);
        
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
    cur_std_withhpc  = std(all_pairs_withhpc(~isnan(all_pairs_withhpc)));
    
    cur_mean_withpfc = mean(all_pairs_withpfc(~isnan(all_pairs_withpfc)));
    cur_std_withpfc  = std(all_pairs_withpfc(~isnan(all_pairs_withpfc)));
    
    mean_corrwithhpc = [mean_corrwithhpc, cur_mean_withhpc];
    mean_corrwithpfc = [mean_corrwithpfc, cur_mean_withpfc];
    
    std_corrwithhpc = [std_corrwithhpc, cur_mean_withhpc];
    std_corrwithpfc = [std_corrwithpfc, cur_mean_withpfc];
end

% Output documentation
% .withhpc_pairs: cell array of size nPatterns, each cell contains a vector

out.withhpc_pairs = withhpc_pairs; % literally list of correlations between source and target=hpc
out.withpfc_pairs = withpfc_pairs; % literally list of correlations between source and target=pfc, but separable for each batch
out.all_pairs_withhpc = all_pairs_withhpc; % same as above, but all in one vector
out.all_pairs_withpfc = all_pairs_withpfc; % same as above, but all in one vector
out.mean_corrwithhpc = mean_corrwithhpc; % mean of all_pairs_withhpc
out.mean_corrwithpfc = mean_corrwithpfc; % mean of all_pairs_withpfc
out.std_corrwithhpc = std_corrwithhpc; % std of all_pairs_withhpc
out.std_corrwithpfc = std_corrwithpfc; % std of all_pairs_withpfc
out.mean_withhpccorr_pattern = mean_withhpccorr_pattern; % mean of withhpc_pairs
out.mean_withpfccorr_pattern = mean_withpfccorr_pattern; % mean of withpfc_pairs
out.prophpc = prophpc; % properties of hpc
out.proppfc = proppfc;

