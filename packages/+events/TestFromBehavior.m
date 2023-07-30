Events = events.initEvents();
% Test case 1: TriggerOn when velocity (vel) exceeds a certain threshold
varnames = {'trajbound'};
triggerOnValues = {1};  % Set the threshold to 2 (outbound)
triggerOnVarnames = {'trajbound'};
Events = events.generateFromBehavior(behavior, varnames, ...
Events, 'triggerOnValues', triggerOnValues, 'triggerOnVarnames', triggerOnVarnames, ...
'windowParameters',[1 1]);
figure;
events.plot_events(Events, 'trajbound');


% Test case 2: TriggerOff when the rat is rewarded (rewarded)
varnames           = {'rewarded'};
triggerOffValues   = {1};  % Set the threshold to 1
triggerOffVarnames = {'rewarded'};
Events = events.generateFromBehavior(behavior, varnames, ...
Events, 'triggerOffValues', triggerOffValues, 'triggerOffVarnames', triggerOffVarnames, ...
'windowParameters',[1 1]);

% Test case 3: TriggerLambda when the rat's velocity (vel) is decreasing (deriv < 0)
varnames          = {'vel', 'deriv'};
triggerVars       = varnames;
triggerVarLambdas = {@(behavior) behavior.deriv < 0, @(behavior) behavior.deriv < 0};
Events = events.generateFromBehavior(behavior, varnames, ...
Events, 'triggerVars', triggerVars, 'triggerVarLambdas', triggerVarLambdas, ...
'windowParameters',[1 1]);

