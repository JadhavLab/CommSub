% Query string
% TABLE : ACQUIRE RUNS
% ----------------------

% Determine keys to use : you can use this string to arbitrarily select rows
%   each item of the filtstring is a property to select. $x pulls the x
%   column and applies the test shown
filtstring = ["ismember($animal, [""JS21"",""ZT2"",""ER1"",""JS14"",""JS13"",""JS17""])",...
    "string($generateH(:,2)) ==  ""fromCoherence  fromRipTimes""",...
    "$spikeBinSize==0.15",...
    "$numPartition==50",...
    "$quantileToMakeWindows == 0.85",...
    "all(cat(1,$winSize{:})==[-0.15,0.15],2)"];

% Get the proper keys
T = table.get.summaryOfRuns('mGH', true);
T = query.table.getHashed_stringFilt(T,filtstring);
[Patterns, otherData] = query.pattern.getHashed_stringFilt(filtstring, 'server', false, 'mGH', true);

[nAnimal, nMethods, nPartition, ~, nResult] = size(Patterns);
disp(['Found ', num2str(nResult), ' results.']);
disp("Results = nAnimal, nMethods, nPartition, nSource, nTarget");
nPatterns = nResult/2;

% Get the number of source and target neurons
nSource = zeros(1,nAnimal);
nTarget = zeros(1,nAnimal);
numDimsUsedForPrediction = cell(1,nAnimal);
for a = 1:nAnimal
    nSource(a) = size(Patterns(a,1,1,1,1).X_source,1);
    nTarget(a) = size(Patterns(a,1,1,1,1).X_target,1);
    numDimsUsedForPrediction{a} = 1:min(nSource(a), nTarget(a));
end

Option = otherData{1}.Option;
if Option.sourceArea == "CA1"
    source = "hpc";
    target = "pfc";
    hpc = 1;
    pfc = 2;
else
    source = "pfc";
    target = "hpc";
    hpc = 2;
    pfc = 1;
end

tempPatterns = permute(Patterns, [2,1,3,4,5]);
newSize = size(tempPatterns);
newSize = [newSize(1), prod(newSize(2:3)), newSize(4:end)];
Patterns_AllAnimals = reshape(tempPatterns, newSize);

nAnimalPartition = nAnimal * nPartition;
patternnames = ["theta", "delta", "ripples"];
directionality = ["hpc-hpc","hpc-pfc"];
patternNames = ["theta", "delta", "ripples","theta-control", "delta-control", "ripples-control"];
genH           = shortcut.generateH(Option.generateH);
direct         = shortcut.directionality(directionality);
patternSymbols = shortcut.patternSymbols(patternNames, 2);
