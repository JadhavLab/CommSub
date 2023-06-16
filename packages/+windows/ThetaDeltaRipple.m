function [cellOfWindows, cutoffs] = ThetaDeltaRipple(Events, Option)

THETA  = 1;
DELTA  = 2;
RIPPLE = 3;

nPatterns = numel(Option.patternNames);

%%--------------------
%% 1s: THETA AND DELTA
%%--------------------
[cellOfWindows, cutoffs] = windows.make(Events.times, Option.quantileToMakeWindows, Events.H(:,THETA:DELTA), Option.winSize);

%%----------------
%% 1b :RIPPLES
%%----------------
if any(contains(Option.generateH, ["fromCoherence","fromWpli"]))
    [cellOfWindows(RIPPLE), cutoffs(RIPPLE)] = windows.make(Events.times,...
        Option.quantileToMakeWindows,  Events.H(:,RIPPLE), Option.winSize,...
        'quantile', Events.Hvals(:,RIPPLE),'higherThanQuantile', true); % % RY: quantile needs to be hvals for ripple coherence/wpli threshold to be correct, but timesd computed from Events.H such that non-ripple times thrown out
else
    [cellOfWindows(RIPPLE), cutoffs(RIPPLE)] = windows.make(Events.times, ...
        1,   Events.H(:,RIPPLE), Option.winSize, 'threshold', 'raw','higherThanQuantile', true);
end

windows.countMessage(cellOfWindows, Option.patternNames,...
    'message', 'initial window creation')

%%----------------------
%% 2: EQUALIZE N OF WINDOWS
%%----------------------
if Option.equalWindowsAcrossPatterns == true
    cellOfWindows = windows.equalizeWindowsAcrossPatterns(cellOfWindows);
end
windows.countMessage(cellOfWindows, Option.patternNames, ...
    'message', 'equalize windows')


numWindowsCut = size(cellOfWindows{1},1);
%windows.printWindowOverlap(cellOfWindows, Option.patternNames);

% -----------------------
% Controls: THETA AND DELTA
% -----------------------
if Option.lowerControl
    quantileControl = 1 - Option.quantileToMakeWindows;
else
    quantileControl = Option.quantileToMakeWindows;
end
if Option.oldControlBehavior
    Hc =  control.generatePatternShuffle(Events.H(:,1:3), Events.times, cellOfWindows); % add control patterns;
else
    Hc = Events.H;
end
[Hc_cellOfWindows, Hc_cutoffs] = windows.make(Events.times,  quantileControl,...  % add windows of control patterns
    Hc(:,THETA:DELTA), Option.winSize,... % Selects less than quantile
    'higherThanQuantile', Option.oldControlBehavior);

% -----------------------
% Controls: RIPPLES
% -----------------------
if any(contains(Option.generateH, ["fromCoherence","fromWpli"]))
    [Hc_cellOfWindows(RIPPLE), Hc_cutoffs(RIPPLE)] = windows.make(Events.times,...
        quantileControl,  Hc(:,RIPPLE), Option.winSize,...
        'quantile', Events.Hvals(:,RIPPLE),'higherThanQuantile', Option.oldControlBehavior);
else
    if Option.oldControlBehavior
        [Hc_cellOfWindows(RIPPLE), Hc_cutoffs(RIPPLE)] = windows.make(Events.times, ...
            1,   Hc(:,RIPPLE), Option.winSize, 'threshold', 'raw','higherThanQuantile', true); %RY quantile won't work because these are raw
    else
        [Hc_cellOfWindows(RIPPLE), Hc_cutoffs(RIPPLE)] = windows.make(Events.times, ...
            quantileControl,   Events.Hvals(:,RIPPLE), Option.winSize, 'threshold', 'quantile','higherThanQuantile', Option.oldControlBehavior); %RY quantile won't work because these are raw
    end
end

% ----------------------------------
% Clean, Remove overalapping windows
% ----------------------------------
% clean up control windows: remove each control pattern's window's overlap
for pattern = 1:length(cellOfWindows)
    curr = windows.removeOverlapsBetweenPattern(...
        cell2mat(cellOfWindows(:,pattern)), cell2mat(Hc_cellOfWindows(:,pattern)));
    Hc_cellOfWindows{pattern} = curr;
end

% -----------------------------------
% Merge many controls into 1 control?
% -----------------------------------
% % Merge into one
cellOfWindows(length(cellOfWindows)+1:length(Hc_cellOfWindows)+length(cellOfWindows)) = Hc_cellOfWindows;
cutoffs = [cutoffs,Hc_cutoffs];

% -----------------------------------------
% Ensure each pattern has equal # of window
% -----------------------------------------
% Equalize trials/windows for each pair of patttern-controlPattern
[cellOfWindows, warnedEmptyControls] =...
    control.equalizePatternControl(cellOfWindows);

% pick control pattern that actually contains controls, would break if all
% three are empty...
if Option.singleControl && warnedEmptyControls
    for iPossibleControl = nPatterns+1:nPatterns*2
        if ~isempty(cellOfWindows{iPossibleControl})
            cellOfWindows{nPatterns+1} = cellOfWindows{iPossibleControl};
            disp("here")
        end
    end
end

