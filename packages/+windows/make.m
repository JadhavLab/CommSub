function [cellOfWindows, cutoff] = make(times, threshold, H, winSize, varargin)
% H: m*n matrix - m is number powers taken (aligned to times, same
%    length as times) and n is the number of rhythms we are interested in
%
% times: time points where the powers of different rhythms were taken
%
% threshold: percentile for which the powers are chosen
%
% winSize: length of window per oscillation
%
% if Opt.quantile given and H, then the Opt.quantile sequencec is used to find
% the quantile() and H is the sequence that the cutoff is applied to in order
% to attain windows.


% Optional argumuents
% -------------------
ip = inputParser;
ip.addParameter('intraWindow_binCenters', []);       % How much to shift between samples .. if not given, it's the winSize
ip.addParameter('intraWindow_binCenter_edges', []);         % How much to shift between samples .. if not given, it's the winSize
ip.addParameter('thresholdMethod','quantile'); % how do we interpret the numbers given in threshold, 'quantile' as a quantilee, and 'raw' as a raw threshold
ip.addParameter('quantile',[]);                % RY if not empty, then it calculates quantile of this variable instead of for H
ip.addParameter('higherThanQuantile', true);   % when there is a single threshold value, this encodes the directionality of the test, >  if true, < if false
ip.parse(varargin{:})
Opt = ip.Results;

boolcheck = [isempty(Opt.intraWindow_binCenters), ...
             isempty(Opt.intraWindow_binCenter_edges)];
assert(all(boolcheck==0) || all(boolcheck==1),...
    'Either provide both intraWindow_binCenters and intraWindow_binCenter_edges or none');

%% 1. calculate the quantiles
%% --------------------------
[nTime,nPatterns] = size(H);
cutoff = zeros(numel(threshold),nPatterns);

% Compute quantile cutoffs either from H or another given variable
for i = 1:nPatterns
     switch Opt.thresholdMethod
        case 'quantile'
            if ~isempty(Opt.quantile)
                cutoff(:,i) = quantile(Opt.quantile(:,i),threshold);
            else
                cutoff(:,i) = quantile(H(:,i),threshold);
            end
        case 'raw'
            cutoff(:,i) = threshold;
        otherwise
            error('Invalid method')
    end
end

%% --------------------------------
%% 2. Apply quantiles to the powers
%% --------------------------------
%marking the H_matrix with 1 if the power is above quantile and 0 if
LOW = 1;
HIGH = 2;
for patternCutoff = 1:length(cutoff)
    for time = 1:nTime
        if size(cutoff, 1) == 2
            if H(time,patternCutoff)>=cutoff(LOW,patternCutoff) & H(time,patternCutoff) < cutoff(HIGH,patternCutoff)
                H(time,patternCutoff)=1;
            else
                H(time,patternCutoff)=0;
            end
        elseif size(cutoff, 1) == 1
            if ~Opt.higherThanQuantile
                if H(time,patternCutoff) < cutoff(patternCutoff)
                    H(time, patternCutoff)=1;
                else
                    H(time, patternCutoff)=0;
                end
            else
                if H(time,patternCutoff) >= cutoff(patternCutoff)
                    H(time, patternCutoff)=1;
                else
                    H(time, patternCutoff)=0;
                end
            end
        else
            error("Cutoff size must be a range (2 numbers) or a min or max (1 number)")
        end
    end
end

if isscalar(winSize)
    winSize = [-winSize, winSize];
elseif iscell(winSize)
    winSize = [winSize{:}];
end


%% -------------------------
%% 3. map to the time points
%% -------------------------
cellOfWindows = cell(1, length(cutoff));
% the time points we want, 1st column starting and 2nd ending
for pattern = 1:length(cutoff)

    switches = diff(H(:,pattern));
    switches = (switches == 1);
    temp = find(switches);
    temp = (switches == 1);
    if iscolumn(times)
        times = times';
    end
    
    R = times(temp)'; % Capture triggers
    starts = R+winSize(1);
    stops  = R+winSize(2);
    if isempty(Opt.intraWindow_binCenters) % Cell of windows is windows x 2 (starts and ends of whole window)
        R = [starts, stops];
    else % Cell of windows is windows x 2 x timeChunks (starts and ends of chunks)
        shifts = 0 : Opt.intraWindow_binCenters : sum(abs(winSize));
        chunkCenters = starts + shifts;
        R = chunkCenters + Opt.intraWindow_binCenter_edges;
    end
    cellOfWindows{pattern} = R;

end
