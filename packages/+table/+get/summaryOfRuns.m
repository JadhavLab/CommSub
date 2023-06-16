function T = summaryOfRuns(varargin)

M = table.matfile.summaryOfRuns(varargin{:});
fields = setdiff(string(fieldnames(M)),"Properties");
if isempty(fields)
    error("Empty data structure!")
end
T = M.(fields(1));
