%function MakeDataFile(varargin)
%
%ip = inputParser;
%ip.addParameter('epoch_type','run');
%ip.parse(varargin{:});
%ip.parse(varargin{:});
%opt = ip.Results;
opt.epoch_type = 'run';
opt.fieldstr   = 'speccoh';

warning off; 
resultsFolder = project('seqNMF');
warning on

%skipMaster = true;

[opt] = seq.initialize(...
               'epoch_type', opt.epoch_type,...
               'fieldstr',   opt.fieldstr,...
               'subsample', 0.085, ... % take a contiguous fraction like this of everyt animal/epoch
               'assignInCaller', false ...
               );

folder = seqnmf_folder(opt, 'usedate', false);
[data, pre] = GetData(opt);
opt.L =  round(timescale/data.params.movingwin(2));

X = data.data;
trainNEURAL = data.data;
testNEURAL  = pre.data(:, ~data.subsample_logical);

mfileData = @(resultsFolder, paramFolder) ...
            matfile(fullfile(resultsFolder, 'datafile_' + string(folder)),...
            'Writable', true);

% Store in mfile
mfile = mfileData(resultsFolder, folder);
mfile.data          = data;
mfile.pre           = pre;
mfile.X             = X;
mfile.trainNEURAL   = trainNEURAL;
mfile.testNEURAL    = testNEURAL;
