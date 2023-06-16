function S = generateRawStruct(varargin)
% Creates the raw structure that summarizes all of the major raw
% data from a given TheScript run.
%
% Standardizes the way raw structures look

ip = inputParser;
ip.KeepUnmatched = true;
ip.addArgument('legacy', false, @islogical);
ip.parse(varargin{:});
Opt = ip.Results;
U = ip.Unmatched;

% Everything that is unmatched seeds the
% raw data structure. But from there, the
% parsing depends on legacy versus non-legacy

if Opt.legacy
    S = struct();
    for field = string(fieldnames(U))
        if startsWith(field, 'X_')
        elseif ismember(field, {'H', 'Hvals', 'Hnanlocs', 'Htimes'})
            S.H.(field) = U.(field);
        else
            S.(field) = U.(field);
        end
    end
else
    for field = string(fieldnames(U))
        switch field
        % Top level
        case {'avgFR', 'areaPerNeuron', 'nPattern'}
            S.(field) = U.(field);
        % .pfc
        case {'nPFCneurons'}
            S.pfc.(field) = U.(field);
        % .hpc
        case {'nHPCneurons'}
            S.hpc.(field) = U.(field);
        % .windowInfo
        case {'cellOfWindows', 'nWindows'}
            S.windowInfo.(field) = U.(field);
        otherwise
            S.(field) = U.(field);
    end

    if all(ismember({'nPatternAndControl', 'nPattern'}, fieldnames(S)))
        S.nControl = S.nPatternAndControl - S.nPattern;
    end
end

