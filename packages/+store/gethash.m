function hash = gethash(Option)
% GETHASH returns a hash of the input string
%
%   Input:
%       Option  - The Option struct assigned for the current run

Option    = nd.flexrmfield(Option, {'loadifexists'});
hash      = DataHash(Option);
hash      = hash(1:7); % Take the first 7 letters of the hash
hash      = string(hash);
