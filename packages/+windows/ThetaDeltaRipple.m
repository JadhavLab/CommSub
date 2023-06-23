function [cellOfWindows, cutoffs] = ThetaDeltaRipple(Events, Option)
% [cellOfWindows, cutoffs] = ThetaDeltaRipple(Events, Option)
%   Generates windows for theta, delta, and ripple events
%   INPUTS:
%       Events: struct with fields:
%           times: 1 x nEvents vector of event times
%           H: nEvents x nPatterns matrix of H values
%           Hvals: nEvents x nPatterns matrix of H values
%       Option: struct with fields:
%           patternNames: 1 x nPatterns cell array of pattern names
%               (e.g. {'SWR', 'Ripple', 'Theta', 'Delta'})
%           winSize: size of window in seconds
%               (e.g. 0.5)
%           quantileToMakeWindows: quantile of H values to make windows
%               (e.g. 0.9)
%           equalWindowsAcrossPatterns: true/false
%               meaning if true, equalizes number of windows across patterns
%           generateH: cell array of strings of H to generate
%               meaning if {'fromCoherence', 'fromWpli'}, then will generate
%               H from coherence and wpli
%           lowerControl: true/false
%               meaning if true, will generate control windows from lower
%               quantile of H values
%           oldControlBehavior: true/false
%               meaning if true, will generate control windows from lower
%               quantile of H values
%           singleControl: true/false
%               meaning if true, will generate control windows from lower
%               quantile of H values
%   OUTPUTS:
%       cellOfWindows: 1 x nPatterns cell array of windows
%       cutoffs: nPatterns x 1 vector of cutoffs

disp("Generating windows for events")
tic;

THETA  = 1;
DELTA  = 2;
RIPPLE = 3;

nPatterns = numel(Option.patternNames);

%%--------------------
%% 1s: THETA AND DELTA
%%--------------------
[cellOfWindows, cutoffs] = windows.make(Events.times, ...
    Option.quantileToMakeWindows, Events.H(:,THETA:DELTA), Option,...
    'positiveDerivativeCheck', true, ...
    'outlierQuantile', Option.thetadelta_outlierQuantile);

%%----------------
%% 1b :RIPPLES
%%----------------
if any(contains(Option.generateH, ["fromCoherence","fromWpli"]))
    [cellOfWindows(RIPPLE), cutoffs(RIPPLE)] = windows.make(Events.times,...
        Option.quantileToMakeWindows,  Events.H(:,RIPPLE), Option,...
        'quantile', Events.Hvals(:,RIPPLE),'higherThanQuantile', true); % % RY: quantile needs to be hvals for ripple coherence/wpli threshold to be correct, but timesd computed from Events.H such that non-ripple times thrown out
else
    [cellOfWindows(RIPPLE), cutoffs(RIPPLE)] = windows.make(Events.times, ...
        1,   Events.H(:,RIPPLE), Option, 'threshold', 'raw','higherThanQuantile', true);
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
[Hc_cellOfWindows, Hc_cutoffs] = ...
    windows.make(Events.times,  quantileControl,...  % add windows of control patterns
    Hc(:,THETA:DELTA), Option,... % Selects less than quantile
    'outlierQuantile', Option.thetadelta_outlierQuantile,...
    'positiveDerivativeCheck', true, ...
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
            1,   Hc(:,RIPPLE), Option, 'threshold', 'raw','higherThanQuantile', true); %RY quantile won't work because these are raw
    else
        [Hc_cellOfWindows(RIPPLE), Hc_cutoffs(RIPPLE)] = windows.make(Events.times, ...
            quantileControl,   Events.Hvals(:,RIPPLE), Option, 'threshold', 'quantile','higherThanQuantile', Option.oldControlBehavior); %RY quantile won't work because these are raw
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

% -------------------------------
% If any cellOfWindows is empty, throw an error
% -------------------------------
if any(cellfun(@isempty, cellOfWindows))
    error("A network pattern has no windows. Check your settings...")
end

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

% -------------------------------
% If any cellOfWindows is empty, throw an error
% -------------------------------
if any(cellfun(@isempty, cellOfWindows))
    error("A network pattern has no windows. Check your settings...")
end

disp("")
disp("Windows generated " + num2str(numWindowsCut) + " windows cut")
disp("  time elapsed: " + num2str(toc) + " seconds")
disp("  cutoffs: " + num2str(cutoffs))
