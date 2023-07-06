function Out = match_rrr(Patterns, Option, Spk, varargin)
% Components = timeVarying(Patterns, Option)
%
% This function takes the patterns and performs the time varying analysis
% on them. It returns a struct with the critical components for each
% behavior and the subspaces for each behavior.
%
% Inputs:
% -----
%   Patterns: struct with the patterns for each partition
%   Option: struct with the options for the analysis
%
% Outputs:
%   Components: struct with the critical components for each behavior and
%   the subspaces for each behavior

ip = inputParser();
ip.parse(varargin{:});
Opt = ip.Results;

disp("Running time varying analysis")
tic
Const = option.constants();

if ~isfield(Patterns, 'rankRegress') || isempty(Patterns(1).rankRegress)
    error("Patterns must have the rankRegress field")
end

for i = progress(1:numel(Patterns), 'Title', 'Time varying analysis - rrr')

    P = Patterns(i); 
    if i > 1
        O = Out(i);
    else
        O = struct();
    end
    O.directionality = P.directionality;
    tmp = split(P.directionality, '-');
    O.source = tmp{1};
    O.target = tmp{2};
    O.match = analysis.temporal.single_match_rrr(P, ...
                Option, Spk);
    O.name = P.name;

    if i == 1
        disp("Initializing output struct")
        Out = repmat(nd.emptyLike(O), size(Patterns));
    end
    Out(i) = O;
end

disp("Finished time varying analysis in " + string(toc) + " seconds")
