function plot(result, varargin)

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('fieldset', {'C','wpli','S1','S2'});
ip.addParameter('ylim', [])
ip.addParameter('smooth', "")
ip.addParameter('interp', "")
ip.addParameter('interpN', [200, 400])
ip.addParameter('qClim',[0.02, 0.98]);
ip.addParameter('lindistlabel',true);
ip.addParameter('cmocean_c', 'balance');
ip.addParameter('cmocean_s', 'curl');
ip.addParameter('unify', false);
ip.addParameter('S1_replace', 'CA1');
ip.addParameter('S2_replace', 'PFC');
ip.addParameter('in',false);
ip.addParameter('yscale','linear');
ip.parse(varargin{:})
args = ip.Unmatched;
opt = ip.Results;

kws = {'FontSize',14};


if contains(opt.smooth,'x')
    for f = opt.fieldset
        try
            result.(f{1}) = smoothdata(result.(f{1}), 'rlowess', 1);
        catch ME
            warning('Cannot smooth field=%s', f{1});
        end
    end
end
if contains(opt.smooth,'y')
    for f = opt.fieldset
        result.(f{1}) = smoothdata(result.(f{1}), 'rlowess', 2);
    end
end
if contains(opt.interp,'x')
    for f = opt.fieldset
        points = 1:size(result.(f{1}),1);
        query_points = linspace(1,size(result.(f{1}),1), opt.interpN(1));
        result.(f{1}) = interp1(points, result.(f{1}), query_points);
    end
end
if contains(opt.interp,'y')
    for f = opt.fieldset
        points        = 1:size(result.(f{1}),2);
        query_points  = linspace(1,size(result.(f{1}),2), opt.interpN(2));
        result.(f{1}) = interp1(points, result.(f{1})', query_points)';
    end
end
if ~isempty(opt.ylim)
    filt = result.f >= opt.ylim(1) & result.f <= opt.ylim(2);
    for f = opt.fieldset
        result.(f{1}) = result.(f{1})(:,filt);
    end
    result.f = result.f(filt);
end

fcnt = 0;
nFields = numel(opt.fieldset);
for field = string(opt.fieldset)
    fcnt = fcnt + 1;
    subplot(nFields,1,fcnt);
    if string(opt.yscale) == "log"
        im=surf(1:size(result.(field),1), linspace(min(result.f), max(result.f), size(result.(field),2)), result.(field)'); 
        set(gca,'zscale', opt.yscale);
        view(2);
        
    else
        im=imagesc(1:100, result.f, result.(field)'); 
    end
    axis xy;
    q = quantile(im.CData(:), opt.qClim);
    set(gca,'clim',q);
    if field.startsWith('S')
        cmocean(opt.cmocean_s)
        if field == "S1"
            field = opt.S1_replace;
        elseif field == "S2"
            field = opt.S2_replace;
        else
            error("You fucked up")
        end
    else
        cmocean(opt.cmocean_c)
    end
    ylabel(string(field) + ' Hz',kws{:})
    xticks([])
    c=colorbar;
    c.Label.String = field;
    climAx = get(gca,'clim');
    climAx(1) = 0;
    set(gca,'Clim',climAx);
end

figaxes=findobj(gcf,'type','axes');
linkaxes(figaxes,'xy')
if opt.in
    set(figaxes,'xdir','reverse')
end

if opt.lindistlabel
    if ~isempty(args)
        args = util.struct2varargin(args);
        linearized.lindistlabel(args{:});
    else
        linearized.lindistlabel();
    end
end



%% Graveyard
%-----------


%if opt.unify
%    q = [min(q(:,1)), min(q(:,2))];
%    set(ax(1),'clim',q);
%    set(ax(2),'clim',q);
%else
%    set(ax(1),'clim',q(1,:));
%    set(ax(2),'clim',q(2,:));
%end
