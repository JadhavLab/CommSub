function [data] = ry_run_seqnmf(data, fields, varargin)

% Optional args
ip = inputParser;
ip.addParameter('unpack', true, @islogical)
ip.addParameter('fieldpack_kws', {} , @iscell);
ip.addParameter('seqnmf_kws', ...
    {'L', 10, 'lambdaL1W', 10, 'lambdaL1W', 10, 'maxiter', 1e3}, @iscell);
ip.parse(varargin{:});
opt = ip.Results;

% Basic seqNMF
% ------------
fieldpack_kws = opt.fieldpack_kws;
% Params
L =  round(10/data.params.movingwin(2));
% Pack fields
data.data = seqnmf_packfields(data, fields, fieldpack_kws{:})';
figure(1001);
%[W, H,cost, loadings, power] = seqNMF_gpu(data.data, ...
%    'L', L, 'lambdaL1W', 10, 'lambdaL1W', 10, 'maxiter', 10e3);
fprintf('SeqNMF args = %s\n', string(opt.seqnmf_kws));
[W, H, cost, loadings, power] = seqNMF_gpu(data.data, opt.seqnmf_kws{:});
data.W = W;
data.H = H;
data.cost = cost;
data.loadings = loadings;
data.power = power;

if opt.unpack
    %data = seqnmf_unpackfields(W, H, data, fields, opt.fieldpack_kws{:});
end
