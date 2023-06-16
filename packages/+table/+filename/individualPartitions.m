function filename = individualPartitions(varargin)
% Returns the table that details each partition run through analyses,
% across all batches of partition-sets. If you run batch of 50
% partitions, these tables absorb everything about the analysis
% of those partitions, continually across instances of running
% TheScript.m
%
% Legacy version used poor naming. Hard to remember. This code
% defaults to writing the moree descriptive table names.

% Parse Inputs
% ------------
ip = inputParser;
ip.addParameter('mGH',    true); %multi-generate Hj
ip.addParameter('legacy', false); 
ip.addParameter('nameOnly', true); 
ip.parse(varargin{:});
Opt = ip.Results;

if Opt.legacy
    if Opt.mGH
        filename = 'Tablerow.mat';
    else
        filename = 'Megatable.mat'; 
    end
else
    if Opt.mGH
        filename = 'table_oneRowPerPartition';
    else
        filename = 'table_oneRowPerPartition'; 
    end
end

if ~Opt.nameOnly
    filename = fullfile(table.filename.folder(), filename);
end

