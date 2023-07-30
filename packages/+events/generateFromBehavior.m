function Events = generateFromBehavior(behavior, varnames, Events, varargin)
    % Create an input parser
    p = inputParser;

    if ~exist('Events', 'var') || isempty(Events)
        % If Events is empty, create a new Events struct
        Events = events.initEvents();
    end
    if isempty(Events.times)
        % If events.times is empty, initialize it
        Events.times = behavior.time;
    end

    % Define default values for optional parameters
    defaultWindowParameters = [0, 0];
    defaultTriggerOnValues = {NaN};
    defaultTriggerOnVarnames = {''};
    defaultTriggerOffValues = {NaN};
    defaultTriggerOffVarnames = {''};
    defaultTriggerVars = {''};
    defaultTriggerVarLambdas = {@(x) false};

    % Add optional parameters to the input parser
    addOptional(p, 'windowParameters', defaultWindowParameters);
    addOptional(p, 'triggerOnValues', defaultTriggerOnValues);
    addOptional(p, 'triggerOnVarnames', defaultTriggerOnVarnames);
    addOptional(p, 'triggerOffValues', defaultTriggerOffValues);
    addOptional(p, 'triggerOffVarnames', defaultTriggerOffVarnames);
    addOptional(p, 'triggerVars', defaultTriggerVars);
    addOptional(p, 'triggerVarLambdas', defaultTriggerVarLambdas);

    % Parse the inputs
    parse(p, varargin{:});

    % For each varname
    for i = 1:length(varnames)
        % Get the variable values
        varvalues = behavior.(varnames{i});

        % Append variable values to H and Hvals fields
        Events.H = [Events.H, varvalues];
        Events.Hvals = [Events.Hvals, varvalues];
    end

    % Get the window parameters
    windowParameters = p.Results.windowParameters;

    % If TriggerOn values and varnames are specified, trigger events based on the detection of the TriggerOn values
    triggerOnVarnames = p.Results.triggerOnVarnames;
    triggerOnValues   = p.Results.triggerOnValues;
    for i = 1:length(triggerOnVarnames)
        if ~isnan(triggerOnValues{i}) && ~isempty(triggerOnVarnames{i})
            disp("Triggering on " + triggerOnVarnames{i} + " = " + triggerOnValues{i});
            field = triggerOnVarnames{i};
            triggerCondition = ... 
            @(behavior) (behavior.(field)(2:end) == triggerOnValues{i}) & ...
            (behavior.(field)(1:end-1) ~= triggerOnValues{i});  % Assuming 'vel' is the triggerOnVarname
            [triggerOnTimes, triggerOnWindows] = getTriggerTimesAndWindows(behavior, triggerCondition, windowParameters);
            Events.cellOfWindows = [Events.cellOfWindows; {triggerOnWindows}];
            % Record the variable name used to generate the event
            Events.cellOfWin_varnames{end+1} = triggerOnVarnames{i};
        end
    end

    % If TriggerOff values and varnames are specified, create windows on ending the specific values
    triggerOffVarnames = p.Results.triggerOffVarnames;
    triggerOffValues = p.Results.triggerOffValues;
    for i = 1:length(triggerOffVarnames)
        if ~isnan(triggerOffValues{i}) && ~isempty(triggerOffVarnames{i})
            % For the TriggerOff situation:
            disp("Triggering off " + triggerOffVarnames{i} + " = " + triggerOffValues{i});
            field = triggerOffVarnames{i};
            triggerCondition = @(behavior) ...
            (behavior.(field)(2:end) ~= triggerOffValues{i}) & (behavior.(field)(1:end-1) == triggerOffValues{i});  % Assuming 'vel' is the triggerOffVarname
            [triggerOffTimes, triggerOffWindows] = getTriggerTimesAndWindows(behavior, triggerCondition, windowParameters);
            Events.cellOfWindows = [Events.cellOfWindows; {triggerOffWindows}];
            Events.cellOfWin_varnames{end+1} = triggerOffVarnames{i};
        end
    end

    % If TriggerVars and TriggerVarLambdas are specified, trigger events based on the lambda functions
    triggerVars = p.Results.triggerVars;
    triggerVarLambdas = p.Results.triggerVarLambdas;
    for i = 1:length(triggerVars)
        if ~isempty(triggerVars{i}) && isa(triggerVarLambdas{i}, 'function_handle')
            disp("Triggering on " + triggerVars{i});
            % For the TriggerLambda situation:
            triggerCondition = triggerVarLambdas{i};  % Assuming triggerVarLambdas{i} is a function handle that takes as input the behavior data and returns a logical array
            [triggerLambdaTimes, triggerLambdaWindows] = getTriggerTimesAndWindows(behavior, triggerCondition, windowParameters);
            Events.cellOfWindows = [Events.cellOfWindows; {triggerLambdaWindows}]; 
            Events.cellOfWin_varnames{end+1} = triggerVars{i};
        end
    end
end

function [triggerTimes, triggerWindows] = getTriggerTimesAndWindows(behavior, triggerCondition, windowParameters)
    % Check the condition for triggering the event
    triggerEvents = triggerCondition(behavior);

    % Get the times when the events are triggered
    triggerTimes = behavior.time([false; triggerEvents]);  % Assuming 'time' is a variable in the behavior table

    % Calculate start and end times for each window
    if windowParameters(1) > 0
        startTimes = triggerTimes - windowParameters(1);
    else
        startTimes = triggerTimes + windowParameters(1);
    end
    endTimes = triggerTimes + windowParameters(2);
    triggerWindows = [startTimes, endTimes];
end

