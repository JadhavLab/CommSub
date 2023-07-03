addpath(hashdefine);

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
Truns = table.get.summaryOfRuns('mGH', true);
Truns = query.table.getHashed_stringFilt(Truns, filtstring);
Tpartitions = table.get.summaryOfPartitions('mGH', true);
Tpartitions = query.table.getHashed_stringFilt(Tpartitions, [filtstring(1), filtstring(3:end-1)]);

keys = unique(Truns.hash);

for key = keys

    load(key + ".mat")
    hash = DataHash(Option);
    hash = hash(1:8); 
    hash = string(hash);
    timestamp = T(hash ==T.hash,:);

    for k = 1:numel(Option.generateH)
        % Get columns for the option struct
        OptionToStore = rmfield(Option, 'generateH');
        OptionToStore.generateH = Option.generateH{k};
        Optiontable = struct2table(OptionToStore, 'AsArray', true);
        
        % Measure how effective this run was
        [optResult, perf] = params.optmizeOptions...
            (theta_sp(k), theta_ripple(k), Patterns(k,:,:,:), nTarget, nSource, Option.numPartition);
        windowCutOff      = cutoffs(k,:);
        theta_sp_corr     = theta_sp(k);
        theta_ripple_corr = theta_ripple(k);
        optResultsrow = table(windowCutOff,theta_sp_corr, theta_ripple_corr ,...
            perf, optResult);
        
        % Measure how effective this run was
        tablerow                   = [table(timestamp, hash), optResultsrow]; % Combine option columnns with hash and date
        currPatterntable           = query.getPatternTable(Patterns(k,:,:,:), Option);
        currPatterntable.generateH = [];
        currPatterntable           = [currPatterntable, repmat(tablerow,height(currPatterntable),1)];
        Patterntable               = [Patterntable; currPatterntable];
    end

    % Which rows to delete?

    % Append

end
