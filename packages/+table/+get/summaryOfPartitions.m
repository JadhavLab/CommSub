function T = summaryOfPartitions(varargin)

M = table.matfile.summaryOfPartitions(varargin{:});
fields = setdiff(string(fieldnames(M)),"Properties");
T = M.(fields(1));

