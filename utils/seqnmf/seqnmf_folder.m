function folder = seqnmf_folder(opt, varargin)

% Parse optional arguments
ip = inputParser();
ip.addParameter('usedate', true)
ip.KeepUnmatched = 1;
ip.parse(varargin{:});
folderopt = ip.Results;

% Setup folder
% ------------
if ~isempty(opt.epoch_type)
    folder = [char(opt.epoch_type) '.'];
else
    folder = '';
end

if ~isempty(opt.fieldstr)
    folder = [folder 'seqnmf_' opt.fieldstr];
    disp(['seqnmf with ' opt.fieldstr]);
end
%if ~isempty(opt.epoch_type) && opt.epoch_type
%    folder = [folder '_' opt.epoch_type];
%end
if ~isempty(opt.timescale) && opt.timescale ~= 0
    opt.lambdafolder = [folder '_' num2str(opt.timescale)];
end
if ~isempty(opt.lambda) && opt.lambda ~= 0
    folder = [folder '_lambda=' sprintf('%1.0e',opt.lambda)];
end
%mark_k = false;
%if opt.timescale > 5
%    if opt.K ~= 4 && ~isequal(fieldstr,'traj')
%        mark_k = true;
%    end
%    if opt.K ~= 8 && isequal(fieldstr,'traj')
%        mark_k = true;
%    end
%else
%    if opt.K ~= 3 && ~isequal(fieldstr,'traj')
%        mark_k = true;
%    end
%    if opt.K ~= 6 && isequal(fieldstr,'traj')
%        mark_k = true;
%    end
%end
if ~isempty(opt.K)
    folder = [folder '_K=' num2str(opt.K)];
end

if ~isempty(opt.lambdaOrthoH) && opt.lambdaOrthoH ~= 0
    folder = [folder '_orthoH=' sprintf('%1.0e', opt.lambdaOrthoH)];
end
if ~isempty(opt.lambdaOrthoW) && opt.lambdaOrthoW ~= 0
    folder = [folder '_orthoW=' sprintf('%1.0e', opt.lambdaOrthoW)];
end

if folderopt.usedate
    folder = [folder '_' date];
end

if isstring(folder)
    folder = join(folder,'');
end
