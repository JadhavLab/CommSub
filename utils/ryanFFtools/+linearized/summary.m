function summary(result, bands, varargin)
% Gramm powered summaries

ip = inputParser;
ip.parse(varargin{:})
opt = ip.Results;

if startsWith(string(opt.cmap), "cmocean_")
    opt.cmap = replace(opt.cmap,'cmocean_', '');
    cmapFunction =  @cmocean;
else
    cmapFunction = @colormap;
end
colors = cmapFunction(opt.cmap, size(bands,1));

if contains(opt.smooth,'x')
    for f = {'C','wpli','S1','S2'}
        result.(f{1}) = smoothdata(result.(f{1}), 'rlowess', 1);
    end
end
if contains(opt.smooth,'y')
    for f = {'C','wpli','S1','S2'}
        result.(f{1}) = smoothdata(result.(f{1}), 'rlowess', 2);
    end
end
if contains(opt.interp,'x')
    for f = {'C','wpli','S1','S2'}
        points = 1:size(result.(f{1}),1);
        query_points = linspace(1,size(result.(f{1}),1), opt.interpN(1));
        result.(f{1}) = interp1(points, result.(f{1}), query_points);
    end
end
if contains(opt.interp,'y')
    for f = {'C','wpli','S1','S2'}
        points        = 1:size(result.(f{1}),2);
        query_points  = linspace(1,size(result.(f{1}),2), opt.interpN(2));
        result.(f{1}) = interp1(points, result.(f{1})', query_points)';
    end
end
if ~isempty(opt.ylim)
    filt = result.f >= opt.ylim(1) & result.f <= opt.ylim(2);
    for f = {'C','wpli','S1','S2'}
        result.(f{1}) = result.(f{1})(:,filt);
    end
    result.f = result.f(filt);
end

subplot(4,1,1)
im=imagesc(1:100, result.f, result.wpli'); axis xy;
q = quantile(im.CData(:), opt.qClim);
set(gca,'clim',q);
cmocean('balance')
ylabel('Hz')
xticks([])
c=colorbar;
c.Label.String = 'wpli';
subplot(4,1,2)
im=imagesc(1:100, result.f, result.C'); axis xy;
q = quantile(im.CData(:), opt.qClim);
set(gca,'clim',q);
cmocean('balance')
ylabel('Hz')
xticks([])
c=colorbar;
c.Label.String = 'C';
subplot(4,1,3)
im=imagesc(1:100, result.f, result.S1'); axis xy;
q = quantile(im.CData(:), opt.qClim);
set(gca,'clim',q);
cmocean('balance')
ylabel('Hz')
xticks([])
c=colorbar;
c.Label.String = 'CA1';
subplot(4,1,4)
im=imagesc(1:100, result.f, result.S2'); axis xy;
Q = quantile(im.CData(:), opt.qClim);
set(gca,'clim',q);
xlabel('Linear Distance')
ylabel('Hz')
cmocean('balance')
c=colorbar;
c.Label.String = 'PFC';


linkaxes(findobj(gcf,'type','axes'),'xy')

