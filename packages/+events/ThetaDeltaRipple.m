function Events = ThetaDeltaRipple(Option)
% Generates the event matrix for theta, delta and ripple events
% 
% USAGE : Events = ThetaDeltaRipple(Option)
%
% INPUT : Option : struct
%   Option.animal : string
%       Name of the animal
%   Option.sourceArea : string
%       Name of the source area
%   Option.generateH : string
%       Method for generating the event matrix
%       "fromSpectra" : Use the spectrogram to generate the event matrix
%       "fromFilteredEEG" : Use the filtered EEG to generate the event matrix
%       "fromCoherence" : Use the coherence to generate the event matrix
%       "fromWpli" : Use the wpli to generate the event matrix
%       "fromRipTimes" : Use the ripple times to generate the event matrix
%   Option.frequenciesPerPattern : int
%       Number of frequencies per pattern
%   Option.patternNames : cell array of strings
%       Names of the patterns
%   Option.globalrippleWindow : int
%       Window size for the global ripple
%   Option.globalrippleWindowUnits : string
%       Units for the window size for the global ripple
%
% OUTPUT : Events : struct
%   Events.times : array of doubles
%       Times of the events, Tx3, where T is the number of events, and 3
%       are the theta, delta and ripple events
%   Events.H : array of doubles
%       Event matrix
%   Events.Hvals : array of doubles
%       Values of the event matrix
%   Events.Hnanlocs : array of doubles
%       Locations of the NaNs in the event matrix

disp("Generating event matrix for " + replace(Option.generateH," ",""));
tic;

const  = Option.shortcut;
RIPPLE = const.RIPPLE;
DELTA  = const.DELTA;
THETA  = const.THETA;

if contains(Option.generateH, "fromSpectra")

    load(Option.animal + "spectralBehavior.mat");
    if Option.sourceArea == "CA1"
        spectrogram   = efizz.S1;
    else
        spectrogram   = efizz.S2;
    end
    
    frequencyAxis = efizz.f;
    times = efizz.t;
    [H, Hvals, Hnanlocs, times] = events.generateFromSpectra(times,...
        spectrogram,...
        frequencyAxis,...
        Option.frequenciesPerPattern);
elseif contains(Option.generateH, "fromFilteredEEG")
    load(Option.animal + "avgeeg.mat");
    [runningSessions, sleepSessions] = getRunningSessions(Option.animal);
    [H, Hvals, Hnanlocs, times] = events.generateFromFilteredEEG(avgeeg, ...
        Option.sourceArea, "patterns",patternNames(1:3),"downsample",10, "sessionsToInclude", runningSessions);
elseif contains(Option.generateH, "fromCoherence")
    load(Option.animal + "spectralBehavior.mat");
    spectrogram = efizz.C;
    frequencyAxis = efizz.f;
    times = efizz.t;
    [H, Hvals, Hnanlocs, times] = events.generateFromSpectra(times, spectrogram, frequencyAxis,...
        Option.frequenciesPerPattern);
elseif contains(Option.generateH, "fromWpli")
    load(Option.animal + "spectralBehavior.mat");
    spectrogram = efizz.wpli;
    frequencyAxis = efizz.f;
    times = efizz.t;
    [H, Hvals, Hnanlocs, times] = events.generateFromSpectra(times, spectrogram, frequencyAxis,...
        Option.frequenciesPerPattern);
    
else
    error("Core method for deriving the event matrix is not recognized");
end

%                                
% ,---.o          |              
% |---'.,---.,---.|    ,---.,---.
% |  \ ||   ||   ||    |---'`---.
% `   ``|---'|---'`---'`---'`---'
%       |    |                   
%% Modify ripple pattern? Lower the threshold as with window sizes
if contains(Option.generateH, "fromRipTimes")

    load(Option.animal + "globalripple01.mat");
    
    if any(contains(Option.generateH, ["fromWpli", "fromCoherence"]))
        [~, H(:,RIPPLE), Hnanlocs(:,RIPPLE), Hvals(:,RIPPLE), minRippleThreshold, ~] = ...
            events.generateFromRipples(globalripple, ...
            'amplitude_at_riptime', true,...
            'rippleBand', Hvals(:,RIPPLE),...
            'rippleBandTime', times,...
            'globalrippleWindowUnits', 'amp');
    else
        % RY: Hvals, not H here for obvious reasons: you  want the original ripple band activity
        [~, H(:,RIPPLE), Hnanlocs(:,RIPPLE), Hvals(:,RIPPLE), minRippleThreshold, ~] = ...
            events.generateFromRipples(globalripple, ...
            'amplitude_at_riptime', true,...
            'rippleBand', Hvals(:,RIPPLE),... 
            'rippleBandTime', times,...
            'globalrippleWindowUnits', 'std');
    end
end

% Assign output struct
Events.times = times;
Events.H = H;
Events.Hvals = Hvals;
Events.Hnanlocs = Hnanlocs;
Events.minRippleThreshold = minRippleThreshold;
disp("Event matrix generated " + toc + " seconds");
