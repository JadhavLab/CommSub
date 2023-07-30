function Events = initEvents(Events)
    if ~exist('Events', 'var') || isempty(Events)
        Events = struct;
        Events.times = [];
        Events.H = []; % the H values are the time series used to create the windows, with nan'd elements for unusable time points
        Events.Hvals = []; % values without unusuable poitns
        Events.Hnanlocs = []; % mask of nan locations nan = unusable, 1 = usable
        Events.minRippleThreshold = [];
        Events.cellOfWindows = {};
        Events.cellOfWin_varnames = {}; % this field represents the variable names that comprise or where used to generate the columns of H. if a single var, a string, else a string vector for multiple creating one H column
        Events.cutoffs = [];
        Events.nWindows = [];
    end
end

