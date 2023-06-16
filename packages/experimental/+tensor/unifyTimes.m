function [A, B, tC] = unifyTimes(A, tA, B, tB, varargin)

ip = inputParser;
ip.addParameter('upsample',true)
ip.addParameter('timedim', -1);
ip.parse(varargin{:})
Opt = ip.Results;

if Opt.upsample
    modify_bit = @(x) x;
else
    modify_bit = @(x) ~x;
end

[~, nAtimes] = size(tA);
[~, nBtimes] = size(tB);

sz_tB = size(tB);
sz_tA = size(tA);

if modify_bit(nBtimes > nAtimes)
    % Interp times
    inds = interp1(tB(:), 1:numel(tB(:)), tA(:), 'nearest');
    % Nans?
    nan_inds = isnan(inds);
    tA = reshape(tA, sz_tA);
    inds(nan_inds) = [];
    tA(nan_inds)   = [];
    tC = tA;
    A = index_tensor(A, find(~nan_inds), Opt);
    % Interp data
    B = index_tensor(B, inds, Opt);
elseif modify_bit(nAtimes >= nBtimes)
    % Interp times
    inds = interp1(tA(:), 1:numel(tA(:)), tB(:), 'nearest');
    tC = tB;
    % Nans?
    nan_inds = isnan(inds);
    inds(nan_inds) = [];
    tB(nan_inds)   = [];
    B = index_tensor(B, find(~nan_inds), Opt);
    % Interp data
    A = index_tensor(A, inds, Opt);
end

function X = index_tensor(X, inds, Opt)

timedim = Opt.timedim;
sz_X = size(X);
if timedim == -1
    timedim = numel(sz_X);
end
X = tens2mat(X, timedim); % mode-n reduction
X = X(inds, :); % interp
sz_X_new = sz_X; % new size
sz_X_new(timedim) = numel(inds); % new dimension
X = mat2tens(X, sz_X_new, timedim); % matrix -> tensor, mode-n
