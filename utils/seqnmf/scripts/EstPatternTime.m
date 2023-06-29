function [bestPattern, pattern] estbestpatterntime(data, W, H)
% ESTBESTPATTERNTIME uses seqnmf pattern H and W to estimate the most prominent
% pattern and temporal location in that pattern at a time
%
% Caveat: Right now only works if one pattern is found
%
% Inputs
% ------
% W    - pattern matrix, NxKxL
% H    - template presence matrix, KxT
% data - matrix of data, NxT
%
% Output
% ------
% bestPattern - list of which pattern L times best match per T
% pattern - list of the total match to each L pattern at each time T

% Test if the user is running improper data through this
[~, Hlocs] = max(H, [], 1);
Hlocs = num2cell(Hlocs);
assert( isequal(Hlocs{:}) );

% Get the relevant sizes of the matrices
L = size(W,3);
T = size(data,2);
N1 = size(data,1);
N2 = size(W, 1);

assert(N2  ==  N1); % Ensure inner dimensions match

% Test each pattern time for the best match
patternTime = zeros(L,T, 'single');
w = squeeze(W(:,1,:));
for l = 1:L
    matchAtL = sum(conv2(data, circshift(w, -l), 'same'), 1); % sum along feature dimension
    patternTime(l, :) = matchAtL;
end

patternLabels = repmat([1:L]', 1, T);

% Identify the best match via the highest maximal inner product match scores



