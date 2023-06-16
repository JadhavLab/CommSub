function Events = ThetaDeltaRipple(Option)

THETA  = 1;
DELTA  = 2;
RIPPLE = 3;

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
if contains(Option.generateH,"fromRipTimes")

    load(Option.animal + "globalripple01.mat");
    
    if any(contains(Option.generateH, ["fromWpli", "fromCoherence"]))
        [~, H(:,RIPPLE), Hnanlocs(:,RIPPLE), Hvals(:,RIPPLE), minRippleThreshold, original] = ...
            events.generateFromRipples(globalripple, ...
            'amplitude_at_riptime', true,...
            'rippleBand', Hvals(:,RIPPLE),...
            'rippleBandTime', times,...
            'globalrippleWindowUnits', 'amp');
    else
        [~, H(:,RIPPLE), Hnanlocs(:,RIPPLE), Hvals(:,RIPPLE), minRippleThreshold, original] = ...
            events.generateFromRipples(globalripple, ...
            'amplitude_at_riptime', true,...
            'rippleBand', Hvals(:,RIPPLE),... RY: Hvals, not H here for obvious reasons: you  want the original ripple band activity
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
