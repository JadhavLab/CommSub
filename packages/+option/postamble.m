function Option = postable(Option)

%% ===============================================
%% Shortcut/alias variables to improve readability
%% ===============================================
animal = Option.animal;
Option.shortcut.THETA  = 1;
Option.shortcut.DELTA  = 2;
Option.shortcut.RIPPLE = 3;
if Option.sourceArea == "CA1"
    Option.shortcut.HPC = 1;
    Option.shortcut.PFC = 2;
else
    Option.shortcut.PFC = 1;
    Option.shortcut.HPC = 2;
end
Option.patternNames = ["theta", "delta", "ripple"];
winSize = Option.winSize(1);
Option.frequenciesPerPattern = [6 10; 0.5 4; 150 200];
[Option.nPatterns,~] = size(Option.frequenciesPerPattern);
if Option.singleControl == true
    Option.nPatternAndControl = Option.nPatterns+1;
else
    Option.nPatternAndControl = Option.nPatterns*2;
end

if ~any(contains(Option.patternNames,"control"))
    Option.patternNames = [Option.patternNames; Option.patternNames+"-control"]';
    Option.patternNames = Option.patternNames(:)';
end

