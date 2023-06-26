function [Patterns] = ccaAnalysis(Patterns, Option)
% ccaAnalysis - Apply CCA to each pattern struct
%   Patterns = ccaAnalysis(Patterns)
%
% INPUTS:
%   Patterns - [1 x 1 struct] see partitionPatterns.m for fields
%
% OUTPUTS:
%   Patterns - [1 x 1 struct] with additional fields:
%       cca - [1 x 1 struct] see cca.core for fields

% Iterate over each pattern and apply CCA
for n = 1:numel(Patterns)
    p = Patterns(n);

    curr_area1 = p.X_source;
    curr_area2 = p.X_target;
    
    % Assuming cca.core function takes two arguments: X_area1 and X_area2
    result = analysis.cca.core(curr_area1, curr_area2);
    p.cca = result;
    
    % Assign the updated struct back to the original array
    Patterns = nd.set(Patterns, n, p);
end
