function out = getHashed_stringFilt(T, filtstring, varargin)
% out = getHashed_stringFilt(T, filtstring, varargin)
%
% Filters a table T using a string filter. The string filter is a string
% that can be evaluated in the context of the table T. For example, if T
% has a column called 'name', then the string filter could be
% 'name == "John"'.
%
% Inputs:
%   T           Table to filter
%   filtstring  String filter
%   mostrecent  If true, only return the most recent entry for each unique
%               combination of variables for the table header columns
%               given

ip = inputParser;
ip.addParameter('mostrecent', {}, @(x) iscell(x) || iscellstr(x) || isstring(x));
ip.parse(varargin{:});
Opt = ip.Results;

if numel(filtstring) > 1
    filtstring = join(filtstring," & ");
end

filtstring = replace(filtstring, '$','T.');
evalstring = sprintf("T(%s, :)", filtstring);
disp("Filtering with:")
disp(evalstring);
disp('');
out = eval(evalstring);

if ~isempty(Opt.mostrecent)
    % Sort the table by the "timestamp" column in descending order
    out = sortrows(out, 'timestamp', 'descend');
    
    % Get the unique combinations of the columns in mostrecent
    [~, ia, ~] = unique(out(:, Opt.mostrecent), 'rows', 'stable');  % Adding 'stable' to preserve order

    out = out(ia, :);
end


end
