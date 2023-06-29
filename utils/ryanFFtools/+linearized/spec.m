function [ax,a]=spec(result, varargin)

ip = inputParser;
ip.addParameter('xlim', []);
ip.addParameter('Slog', false)
ip.addParameter('Szscore', false)
ip.addParameter('SmultFreq', true)
ip.parse(varargin{:})
opt = ip.Results;
kws={'FontSize',14};

if ~isempty(opt.xlim)
    filt = result.f >= opt.xlim(1) & result.f <= opt.xlim(2);
    for f = {'C','wpli','S1','S2'}
        result.(f{1}) = result.(f{1})(:,filt);
    end
    result.f = result.f(filt);
end

colorsC=cmocean('balance',4);
colorsS=cmocean('curl',4);
epsilon = 0.02;

X=nanmean(result.wpli);
subplot(4,1,1)
a(1) = area(result.f, X);
set(gca,'xlim',opt.xlim);
ylabel('WPLI',kws{:})
ylim([max(min(X),0), max(X)+epsilon])
set(a(end), 'FaceColor', colorsC(2,:), 'EdgeColor', colorsC(3,:),'LineWidth',2);

X = nanmean(result.C);
subplot(4,1,2)
a(2)=area(result.f, X);
set(gca,'xlim',opt.xlim);
ylabel(sprintf('Coherence\n(Chronux)'),kws{:})
ylim([min(X), max(X)+epsilon])
set(a(end), 'FaceColor', colorsC(2,:), 'EdgeColor', colorsC(3,:),'LineWidth',2);

X = result.S1;
if opt.Slog
    X = log10(X);
end
if opt.Szscore
    X = bsxfun(@rdivide, X, std(X));
end
if opt.SmultFreq
    X = bsxfun(@times, X, result.f);
end
X=nanmean(X);
subplot(4,1,3)
a(3)=area(result.f, X);
set(gca,'xlim',opt.xlim);
ylabel(sprintf('CA1\nPower'),kws{:})
ylim([min(X), max(X)])
set(a(end), 'FaceColor', colorsS(3,:));

X = result.S2;
if opt.Slog
    X = log10(X);
end
if opt.Szscore
    X = bsxfun(@rdivide, X, std(X));
end
if opt.SmultFreq
    X = bsxfun(@times, X, result.f);
end
X=nanmean(X);
subplot(4,1,4)
a(4)=area(result.f, X);
set(gca,'xlim',opt.xlim);
ylabel(sprintf('PFC\nPower'),kws{:})
ylim([min(X), max(X)])
set(a(end), 'FaceColor', colorsC(3,:), 'EdgeColor', colorsC(2,:),'LineWidth',2);

xlabel(sprintf('Frequency\n(Hz)'),kws{:});

ax = get(gcf,'Children');