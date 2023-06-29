function [imCell, indplot] = ry_WHPlot(W,H,Data, varargin)  
%
% plots seqNMF factors W and H
% Also plots Data if provided, and reconstruction if data is not provided
% plotAll=1 means plot all data
% plotAll=0 means crop data so it's on the same timescale as W's 
%
% Emily Mackevicius 2/3/2018

%% Partse optional inputs
ip = inputParser;
ip.addParameter('plotAll',1); 
ip.addParameter('ExtraMatrixToPlot', []); 
ip.addParameter('plotWflat', 0)
ip.addParameter('extraOverlay', false);
ip.addParameter('extraEnlargeTop', false);
ip.addParameter('defaultPalette', 'lfp')
ip.addParameter('colorschemeW', []);
ip.addParameter('colorscheme',  []);
ip.addParameter('pctscaleExtraMatrix', false);
ip.addParameter('optomizetime', true);
ip.addParameter('imgaussfiltscale', true);
ip.addParameter('static', false)
% --------- Axis options -------------------
ip.addParameter('Fs', 0);             % If present, add a time axis
ip.addParameter('t', []);             % If present, add a time axis
ip.addParameter('xaxis_pts', 5);      % If present, add a time axis
ip.addParameter('yaxis_centers', []); % Pass a y-axis if you want one
ip.addParameter('yaxis_axP', []);     % Pass a y-axis if you want one
% --------- Subplot size options -------------------
ip.addParameter('wwidth',0.2)
% --------- Data rep options -----------
ip.addParameter('Hcenter', true)
ip.addParameter('colorWgroups', {})
ip.addParameter('colorWmaps', 'curl', @(x) iscell(x) || ischar(x))
ip.addParameter('meanWgroup', {})
ip.addParameter('smoothWgroup', {})
ip.addParameter('removeWgroup', {})
% --------- Index options --------------
ip.addParameter('i', []);             % Where to start the plot
ip.addParameter('indplot', []);             % Where to start the plot
ip.parse(varargin{:})
opt = ip.Results;

% Ratios
% ------
A = 2.5; % Normal Buffer zone added for lines drawn
B = 1.6; % Factor of how much to extend for extra data
C = 4;   % Buffer zone added for lines drawn around data
dataplotsratio = A/C;

clf
if nargin < 3
    Data = [];
end

if ~isempty(opt.removeWgroup)
    group = cat(1,opt.removeWgroup{:});
    W(group,:,:) = [];
    Data(group,:) = [];
end

%% Options
% Size!
[N, K, L] = size(W);
[~, T]    = size(H);

[N,K,L] = size(W); 
[~,T] = size(H);
color_palet = [[0 .6 .3]; [.7 0 .7]; [1 .6 0];  [.1 .3 .9];  [1 .1 .1];  [0 .9 .3]; [.4 .2 .7]; [.7 .2 .1]; [.1 .8 1 ]; [1 .3 .7]; [.2 .8 .2]; [.7 .4 1]; [.9 .6 .4]; [0 .6 1]; [1 .1 .3]]; 
switch opt.defaultPalette
case 'spikes'
    if isempty(opt.colorscheme)
        opt.colorscheme = '-grays';
    end
    if isempty(opt.colorschemeW)
        opt.colorscheme = '-grays';
    end
case 'lfp'
%% Options
    set(gcf, 'color', 'w');
% Size!
    color_palet = repmat(color_palet, ceil(K/size(color_palet,1)),1); 
    kColors = color_palet(1:K,:); 
    if isempty(opt.colorscheme)
        opt.colorscheme = color_palet;
    end
    if isempty(opt.colorschemeW)
        opt.colorscheme = color_palet;
    end
end
color_extraMatrix = [158,1,66;213,62,79;244,109,67;253,174,97;254,224,139;255,255,191;230,245,152;171,221,164;102,194,165;50,136,189;94,79,162;1 1 1];
% cleanerW = W.*(W>max(W(:))*.1);
% cleanerH = H.*(H>max(H(:))*.1); 
% set widths of subplots
margin = .05; % margin
wwidth = opt.wwidth; % width of W plot
if opt.plotWflat
    wwflat = .025; % width of Wflat plot
else 
    wwflat = 0;
end
hheight = .20; % height of H plot
hdata = 1-hheight; %1-hheight+1*margin; ;/
wdata = 1-wwidth-wwflat-2*margin; 
sep = ceil(L*.5); 

%% Crop data, unless opt.plotAll
%% -----------
if opt.plotAll
    indplot = 1:T;
elseif ~isempty(opt.indplot)
    indplot = opt.indplot;
    suptitle(sprintf('start=%d, end=%d',opt.indplot(1), opt.indplot(end)));
elseif ~isempty(opt.i)
    indplot = 2*L+(1:ceil((K*(L+sep))/wwidth*wdata)); % so that data and W's are on same scale
    indplot = opt.i:opt.i+indplot;
    suptitle(sprintf('i=%d',opt.i));
else
    indplot = 2*L+(1:ceil((K*(L+sep))/wwidth*wdata)); % so that data and W's are on same scale
    if opt.optomizetime == true && ~isempty(opt.ExtraMatrixToPlot)
        % Find portion of the pattern where the convolved output best matches the raw data
        [~,start] = min( sum((opt.ExtraMatrixToPlot - imgaussfilt(Data)).^2,1) )
    elseif opt.static == false
        start = randi(T-2*L-1,1);
    end
    indplot = indplot+start;
    indplot(indplot>T) = [];
    TE = size(opt.ExtraMatrixToPlot,2); 
    if TE>0
        opt.ExtraMatrixToPlot(:,[1:floor(indplot(1)/T*TE) ceil(indplot(end)/T*TE):end])=[];
    end
end

%% --------------------PLOT DATA--------------------
axIm = subplot('Position', [margin+wwidth margin wdata hdata]);
cla
if ~isempty(Data)
    clims = [0 prctile(Data(Data>0),99)]; 
    %helper.grayscalePatchPlot(flipud(Data(:,indplot)));
    if opt.imgaussfiltscale
        %Im = imagesc(smoothdata(Data(:,indplot), 2, 'SmoothingFactor', 0.5), clims);
       D = imgaussfilt(Data(:,indplot), 1, 'FilterSize', [7 3]);
       D = smoothdata(D);
       Im = imagesc(D, clims);
    else
        Im = imagesc(Data(:,indplot), clims);
    end
    cmocean(opt.colorscheme)
else
    toplot = helper.reconstruct(W,H(:,indplot));
    clims = [0 prctile(toplot(:),99.9)]; 
    cla
    % for ki = 1:K
    % toplot = helper.reconstruct(W(:,ki,:),H(ki,indplot));
    % helper.grayscalePatchPlot(flipud(toplot), [], kColors(ki,:));
    % end
    %helper.grayscalePatchPlot(flipud(toplot));
    Im = imagesc(flipud(toplot),clims);
    cmocean(opt.colorscheme)
end
set(gca,'ydir','normal')
hold on; 
plot([1 1 length(indplot) length(indplot) 1], [0 N N 0 0]+.5, 'k');
xlim([0 length(indplot)+1]); ylim([0 N+A])
if opt.Fs || ~isempty(opt.t)
    ax = gca;
    set(gca,'TickLength',[0.001, 0.001])
    xticks(linspace(ax.XLim(1), ax.XLim(2), opt.xaxis_pts));
    if opt.Fs
        xlabs = 1:numel(indplot);
        xlabs = xlabs * (1/opt.Fs);
        xticklabels(linspace(xlabs(1), xlabs(end), opt.xaxis_pts));
    else
        xlabs = t(indplot);
        xticklabels(linspace(xlabs(1), xlabs(2), opt.xaxis_pts));
    end
end
%yticks( opt.yaxis_centers );
%yticklabels( opt.yaxis_axP );
yticks([]);
if ~isempty(opt.yaxis_centers) && ~opt.Fs && ~isempty(opt.t)
    axis off
end

%% --------------------PLOT EXTRADATA--------------------
if ~isempty(opt.ExtraMatrixToPlot) % there is an extra matrix
    if opt.pctscaleExtraMatrix % scale by percentile?
        opt.ExtraMatrixToPlot = opt.ExtraMatrixToPlot - prctile(opt.ExtraMatrixToPlot(:),50); 
        opt.ExtraMatrixToPlot(opt.ExtraMatrixToPlot<=0)=0;
    end
    if opt.extraOverlay % use an OVERLAY strategy (ie plot on top)
        % Construct an Rgb (red matrix for signal)
        tmp = zeros([size(opt.ExtraMatrixToPlot), 3]);
        tmp(:,:,1) = flipud(opt.ExtraMatrixToPlot(:,indplot));
        tmp = 1 * (tmp-min(tmp(:)))/range(tmp(:));
        alphadata = tmp(:,:,1) * -1 + 1; % AlphaData should be inversion, big signal == small alpha value
        im=image(tmp, 'AlphaData', alphadata);
        im.Tag = 'ExtraData';
    elseif opt.extraEnlargeTop
        warning('Not yet implemented')
    else % Simply append to the empty space above the data
        %im = imagesc(opt.ExtraMatrixToPlot, 'ydata', [N+A B*N+A], 'xdata', [1 length(indplot)]); 
        im = imagesc(opt.ExtraMatrixToPlot, 'ydata', [N+A B*N+A], 'xdata', [1 length(indplot)]); 
        cmap  = 1/256*flipud(color_extraMatrix); % from colorbrewer spectral
        cmap(1,:) = 0; 
        colormap(cmap)
        ylim([0 B*N+C])
        cmocean(opt.colorscheme)

    end
end

%% -------------------- PLOT W'S --------------------
axW = subplot('Position', [margin margin wwidth hdata]);
cla
hold on
set(gca, 'ColorOrder', kColors); 
if opt.extraEnlargeTop % WHETHER TO MAKE A SECOND AXIS TO SHOW WS BY THE EXTRA DATA
    nest = 1;
    i = 2;
else
    i = 1;
end
WsToPlot = zeros(N,K*(L+sep)); 
XsToPlot = zeros(3,K); 
YsToPlot = [zeros(1,K); N*ones(2,K)]+.5;
for ki = 1:K
    XsToPlot(:,ki) = [(L+sep)*(ki-1) (L+sep)*(ki-1) (L+sep)*(ki-1)+L+1];
    component = squeeze(W(:,ki,:));
    if ~isempty(opt.smoothWgroup)
        groupCnt = 0;
        for group = opt.smoothWgroup
            groupCnt = groupCnt + 1;
            group = group{1};
            component(group, :) = smoothdata(component(group,:), 2, 'movmean', round(L/20));
        end
    end
    if ~isempty(opt.meanWgroup)
        groupCnt = 0;
        for group = opt.meanWgroup
            groupCnt = groupCnt + 1;
            group = group{1};
            component(group, :) = repmat(2*mean(component(group,:), 2), 1, size(component(group,:),2));
        end
    end
    WsToPlot(:,((L+sep)*(ki-1)+1):((L+sep)*(ki-1)+L)) = component;
end
clims = [0 prctile(WsToPlot(WsToPlot>0),99)]; % if all W's are empty this line will bug
%imagesc(flipud(WsToPlot), clims);%, clims(2)); 
if ~isempty(opt.colorWgroups)
    % Apply percentile manually
    WsToPlot(WsToPlot> clims(2)) = clims(2);
    % Set default color to grayscale for noncolored groups
    ncolors = 10000;
    colorscale = cmocean(opt.colorschemeW, ncolors);
    originalSize = size(WsToPlot);
    Worig = WsToPlot;
    WsToPlot = discretize(WsToPlot(:), ncolors);
    WsToPlot = colorscale(WsToPlot,:);
    WsToPlot = reshape(WsToPlot, [originalSize, 3]);
    % Now per color group, set a different color per the number of rows
    groupCnt = 0;
    for group = opt.colorWgroups
        group = group{1};
        groupCnt = groupCnt + 1;
        if iscell(opt.colorWmaps)
            tmpcolormap = opt.colorWmaps{groupCnt};
        else
            tmpcolormap = opt.colorWmaps;
        end
        dGroup = diff(group);
        if median(dGroup) > 1
            group = min(group):max(group)+dGroup-1;
        end
        colorscale = cmocean(tmpcolormap, numel(group));
        group = group(ismember(group,1:size(WsToPlot,1))); % Throw out labels that are too far
        Wslice = WsToPlot(group,:,:);
        Worigslice = Worig(group,:);
        for row = 1:size(Worigslice,1)
            for col = 1:size(Worigslice,2)
                if ~Worigslice(row,col) == 0
                    Wslice(row,col,:) = Worigslice(row,col) .* permute(colorscale(row,:),[1 3 2]);
                end
            end
            Wslice(row,:,:) = bsxfun(@rdivide, Wslice(row,:,:), max(Wslice(row,:,:),2));
            for col = 1:size(Worigslice,2)
                if Worigslice(row,col) == 0
                    Wslice(row,col,:) = [1, 1, 1];
                end
            end
        end
        WsToPlot(group,:,:) = Wslice;
        %colorscale(1,:) = [1, 1, 1];
        %Wslice = Worig(group, :);
        %originalSize = size(Wslice);
        %Wslice = discretize(Wslice(:), numel(group));
        %Wslice = colorscale(Wslice, :);
        %Wslice = reshape(Wslice, [originalSize, 3]);
        %WsToPlot(group, :, :) = Wslice;
    end
    image((WsToPlot));%, clims(2)); 
else
    WsToPlot = imgaussfilt(WsToPlot, 1, 'FilterSize', [3 3]);
    imagesc((WsToPlot), clims);%, clims(2)); 
    cmocean(opt.colorschemeW);
end
plot(XsToPlot,YsToPlot, 'linewidth', 4);

% Establish x and y limits
xlim([0 K*(L+sep)]);
if length(opt.ExtraMatrixToPlot)>0
    ylim([0 B*N+C])
    %subplot('Position', [margin (1/B)*hdata wwidth hdata]) % Cover up the extra axis protruding beyond its rightful limit
else
    ylim([0 N+A])
end
if ~isempty(opt.yaxis_centers)
    newMax = N;
    oldMax = size(Data,1);
    opt.yaxis_centers = opt.yaxis_centers * newMax/oldMax;
    yticks(opt.yaxis_centers);
    yticklabels(opt.yaxis_axP);
end
set(gca,'TickLength',[0.001, 0.001])
%set(gca,'YColor',[1 1 1])

set(gca, 'ydir', 'normal')
%axis off

%% ---------- PLOT WFLAT (collapse out L dimension of W)----------
if opt.plotWflat
    axWflat = subplot('Position', [margin+wwidth+wdata margin wwflat hdata]);cla
    hold on
    set(gca, 'ColorOrder', kColors); 
    plot(squeeze(sum(W,3)), 1:1:N,':.', 'markersize', 1);
    %keyboard
    %for row = 1:size(kColors)
    %    Wslice = squeeze(sum(W,3));
    %    area(Wslice(:,row), 1:1:N, 'facealpha', 0.1);
    %end
    axis tight
    if ~isempty(opt.ExtraMatrixToPlot)
        ylim([~ise0 B*N+C])
    else
        ylim([0 N+A])
    end

    xlims = xlim; 
    xlim([xlims(2)*.1 xlims(2)])
    set(gca, 'ydir', 'normal')
    axis off
end

%% --------------------PLOT H's--------------------
axH = subplot('Position', [margin+wwidth margin+hdata wdata hheight]);
cla
Hrescaled = repmat(squeeze(sum(sum(W,1),3))',1,T).*H; % rescale by approximate loading
dn = prctile(Hrescaled(:)+eps,100)/2;
if opt.Hcenter
    shift = L / 2;
else
    shift = 0;
end
for ki = K:-1:1
    Xs = [1 1:length(indplot) length(indplot)]; 
    Ys = [dn*ki (dn*ki + Hrescaled(K-ki+1, indplot-shift)) dn*ki]-dn/2;
    patch(Xs,Ys, kColors(K-ki+1,:), 'edgecolor', kColors(K-ki+1,:))
    hold on
end
ylim([dn/2 dn*K+dn*3]);xlim([0 length(indplot)+1])
axis off

%% ---------- Axis linkage and Drawing ----------
if opt.plotAll
%       linkaxes([axIm axW axWflat], 'y'); 
      linkaxes([axIm axH], 'x');
end
if opt.plotWflat
    linkaxes([axW axIm, axWflat], 'y');
else
    linkaxes([axW axIm], 'y');
end
linkaxes([axIm axH], 'x');
imCell = {axIm, axW, axH};

drawnow

