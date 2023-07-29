function Events = initEvents(Events)
    if ~exist('Events', 'var') || isempty(Events)
        Events = struct;
        Events.times = [];
        Events.H = [];
        Events.Hvals = [];
        Events.Hnanlocs = [];
        Events.minRippleThreshold = [];
        Events.cellOfWindows = {};
        Events.cellOfWin_varnames = {};
        Events.cutoffs = [];
        Events.nWindows = [];
    end
end

