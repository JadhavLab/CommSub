function spectrum_1(data, varargin)
% SPECTRUM_1 Spectrum without Time |_K


[P, T] = size(data);

f = gcf;

plot(1:P, mean(data.data,2));
set(gca, 'ydir','normal');
cmocean(opt.colorscheme);
xticks( axisCenters )
xticklabels( axP );
