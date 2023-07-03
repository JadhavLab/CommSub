function savevars(Option, Events, Spk, varargin)


savestruct = struct('Option', Option, varargin{:});

thisFile = fullfile(hashdefine(), hash + ".mat");
disp("Saving ...");
tic; save(thisFile, saveVars{:},'-v7.3');
disp("... " + toc + " seconds");
% link most recent state
pushd(hashdefine());
recencyName = Option.animal + "_" + replace(Option.generateH," ", "") + ...
	    "_mostRecentState.mat";
system("ln -sf " + hash + ".mat " + recencyName);
popd()
% save raw?
if Option.saveRaw
disp("Saving raw...");
tic; save(thisFile, "Events", "Spk",'-v7.3', '-append');
disp("... " + toc + " seconds");
end
