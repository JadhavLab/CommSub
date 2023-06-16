function [cellOfWindows, cutoff] = make_withPhase(times, threshold, H, winSize, varargin)
%function [cellOfWindows, cutoff] = make_withPhase(times, threshold, H, winSize, varargin)
% H: m*n*2 matrix - m is number powers taken (aligned to times, same
%    length as times) and n is the number of rhythms we are interested in
%    the last dimension comprises amplitude and phase information respectively
%
%    Algorithm searches for the first instance of a target phase inside the
%    window of possible amplitudes, and then uses winSize (in cycles) to
%    determine the start and stop points, shifting by winShiftSize (in cycles)
%
% times: time points where the powers of different rhythms were taken
%
% threshold: percentile for which the powers are chosen
%
% winSize: length of window per oscillation
%
% %%% OPTIONALS %%%%%
% if Opt.quantile given and H, then the Opt.quantile sequencec is used to find
% the quantile() and H is the sequence that the cutoff is applied to in order
% to attain windows.


% Optional argumuents
% -------------------
ip = inputParser;
ip.addParameter('winShiftSize', []);  % How much to shift between samples .. if not given, it's the winSize
ip.addParameter('thresholdMethod','quantile'); % how do we interpret the numbers given in threshold, 'quantile' as a quantilee, and 'raw' as a raw threshold
ip.addParameter('quantile',[]);                % RY if not empty, then it calculates quantile of this variable instead of for H
ip.addParameter('higherThanQuantile', true);   % when there is a single threshold value, this encodes the directionality of the test, >  if true, < if false
ip.parse(varargin{:})
opt = ip.Results;

%% 1. calculate the quantiles
%% --------------------------
[nTime,nPatterns] = size(H);
cutoff = zeros(numel(threshold),nPatterns);

% Compute quantile cutoffs either from H or another given variable
for i = 1:nPatterns
     switch opt.thresholdMethod
        case 'quantile'
            if ~isempty(opt.quantile)
                cutoff(:,i) = quantile(opt.quantile(:,i),threshold);
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
            if ~opt.higherThanQuantile
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

%% -------------------------
%% 3. map to the time points
%% -------------------------
cellOfWindows = cell(1,length(cutoff));
% the time points we want, 1st column starting and 2nd ending
for pattern = 1:length(cutoff)

    switches = diff(H(:,pattern));
    switches = (switches == 1);
    temp = find(switches);
    temp = (switches == 1);
    if iscolumn(times)
        times = times';
    end
    
    R = times(temp)';
    if isscalar(winSize)
        R(:,2) = R(:,1)+winSize;
    elseif numel(winSize) == 2
        R = [R+winSize(1), R+winSize(2)];
    else
        error("winSize either a single number (total window time from the trigger) or two numbers (time in front of and behind the trigger)")
    end
    
    cellOfWindows{pattern} = R;

end
