function prettyLabels = prettyFreqLabel(rawLabels)
% PRETTYFREQLABEL takes numerical frequency and generates pretty labels
prettyLabels = string(arrayfun(@(x) sprintf(['%2.' num2str(round_precision) 'f hz'], x), rawLabels, 'UniformOutput', false));

