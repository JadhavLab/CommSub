
% Find parameters
% ---------------
%seqnmf_findLambda 


% Major analyses
% --------------
nmfs = {'overall'}; % modes : animal_overall (animal-wise overall patterns), overall (patterns for entire set of animals at once)
seqnmf % Run overall
nmfs = {'animal_overall'}; % modes : animal_overall (animal-wise overall patterns), overall (patterns for entire set of animals at once)
seqnmf % Run for animals

% Cross-validation
% ----------------
nmfs = {'overall'}; % modes : animal_overall (animal-wise overall patterns), overall (patterns for entire set of animals at once)
seqnmf_crossval

