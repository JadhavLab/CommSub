
function countMessage(cellOfWindows, patternNames, varargin)

ip = inputParser;
ip.addParameter('message', 'windows');
ip.parse(varargin{:})
Opt = ip.Results;

disp(newline);
disp('--------------------------')
disp(Opt.message)
disp('--------------------------')
cellfun(@(x,y)  fprintf("%d windows for %s\n", size(x,1),  y), ...
    cellOfWindows, cellstr(patternNames(1:3)));
cellfun(@(x,y)  fprintf("%d timerange for %s\n", range(x(:)),  y), ...
    cellOfWindows, cellstr(patternNames(1:3)));
