
% --------------------
%  ____  _   _ _   _ 
% |  _ \| | | | \ | |
% | |_) | | | |  \| |
% |  _ <| |_| | |\  |
% |_| \_\\___/|_| \_|
%                    
% --------------------
tic; fig_name = char("Live " + seqStyle + ": " + ...
                replace(folder, 'seqnmf_',''));
fig(fig_name)
reset(gpuDevice())
L =  round(timescale/data.params.movingwin(2))
%[W, H,cost, loadings, power] = seqNMF_gpu(data.data, ...
%    'L', L, 'lambdaL1W', 10, 'lambdaL1W', 10, 'maxiter', 10e3);
if timescale > 5 
    reconstruct_gpu = 0;
else
    reconstruct_gpu = 1;
end
{fieldstr, timescale, epoch_type}
[W, H, cost, loadings, power] = seqNMF_gpu(data.data, kws{:}, 'gpuReconstruct', reconstruct_gpu, 'L', L, 'fig_name', fig_name);
% Record parameters and results into data
data.W             = W;
data.H             = H;
data.WH = helper.reconstruct(data.W, data.H)
data.L             = L;
data.K             = K
data.cost          = cost;
data.loadings      = loadings;
data.power         = power;
data.seqnmf_kws    = seqnmf_kws;
data.kws           = kws;
data.samprate      = 1/median(diff(data.t));
data.timescale     = timescale;
data.seqStyle      = seqStyle
data.fieldstr      = fieldstr
data.fieldpack_kws = fieldpack_kws;
data = seq.unpackfields(data, F);toc
% Save simpleplot computed during run
simple_plot = fig(fig_name);
SimpleWHPlot(W, H, data.data, 0);
saveas(simple_plot, sprintf('%s.%s', file, 'svg'))
saveas(simple_plot, sprintf('%s.%s', file, 'fig'))

