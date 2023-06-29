function eegdat = applyFilterFiles(eegdat, varargin)
% Used to apply filters to the various eeg files

ip = inputParser;
ip.addParameter('filternames',...
                ["delta","theta","ripple"]);
ip.addParameter('filterfiles',...
                ["deltafilter.mat","thetafilter.mat","ripplefilter.mat"]);
ip.addParameter('filterpath',...
                "~/Code/Src_Matlab/usrlocal/filtering");
ip.parse(varargin{:});
opt = ip.Results;

addpath(genpath(opt.filterpath));
onCleanup(@() rmpath(genpath(opt.filterpath)));

% Get the filters
for f = 1:numel(opt.filternames)
    tmp = load(opt.filterfiles(f));
    filterdata{f} = tmp.(replace(opt.filterfiles(f),'.mat',''));
end

% What indices to apply to?
indices = ndBranch.indicesMatrixForm(eegdat);
eegdat = ndBranch.toNd(eegdat);

for i = 1:numel(eegdat)
    eegdat(i).data = squeeze(eegdat(i).data);
    eegdat(i).data = eegdat(i).data(:,1);
end

% Iterate and apply
tmp = cell(1, numel(opt.filterfiles));
for index = progress(indices','Title','Epochs')

    I = index;
    I = num2cell(I);
    for f = 1:numel(opt.filternames)
        name = opt.filternames(f);
        eegdat( I{:} ).(name) = ...
            filtereeg2(eegdat( I{:} ), ...
                       filterdata{f}, 'int16', 1);
    end

end

% Restore ND shape to BranchedND
eegdat = nd.branched(eegdat);
