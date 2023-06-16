function windows(cellOfWindows, patternNames, message)
disp(newline);
disp('--------------------------')
fprintf('%s:', message);
disp('--------------------------')
if numel(patternNames) > numel(cellOfWindows)
    patternNames = patternNames(1:numel(cellOfWindows));
if numel(patternNames) < numel(cellOfWindows)
    patternNames = patternNames(1:numel(cellOfWindows));
end
N = numel(cellOfWindows);
cellfun(@(x,y)  fprintf("%d windows for %s\n", size(x,1),  y), ...
    cellOfWindows, cellstr(patternNames));
cellfun(@(x,y)  fprintf("%d timerange for %s\n", range(x(:)),  y), ...
    cellOfWindows, cellstr(patternNames));
