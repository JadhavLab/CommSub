function colornames_view(palette,order)
% View the COLORNAMES palettes in an interactive figure. Sort colors by name/colorspace.
%
% (c) 2014-2022 Stephen Cobeldick
%
%%% Syntax:
%  colornames_view
%  colornames_view(palette)
%  colornames_view(palette,order)
%
% Create a figure displaying all of the colors from any palette supported
% by the function COLORNAMES. The palette and sort order can be selected
% by drop-down menu or by optional inputs. The colors may be sorted:
% * alphanumerically (names include any leading indices), or
% * alphabetically (names exclude any leading indices), or
% * by colorspace: Lab, LCh, XYZ, YUV, HSV, or RGB.
%
% Note0: Requires the function COLORNAMES and its associated MAT file (FEX 48155).
% Note1: (re-)drawing the colornames is slow for larger color palettes.
%
%% Input and Output Arguments %%
%
%%% Inputs (all inputs are optional):
%  palette   = CharRowVector, the name of a palette supported by COLORNAMES.
%  sortorder = CharRowVector, the colorspace dimensions in the desired order,
%              eg: 'Lab', 'abL', 'bLa', etc. OR 'Alphabet' OR 'AlphaNum'.
%
%%% Outputs:
% none
%
% See also COLORNAMES COLORNAMES_CUBE COLORNAMES_DELTAE MAXDISTCOLOR COLORMAP

%% Figure Parameters %%
%
persistent fgh axh slh pah srh txh txs cnc rgb coh prv
%
isChRo = @(s) ischar(s) && ndims(s)==2 && size(s,1)==1; %#ok<ISMAT>
%
% Text lightness threshold:
thr = 0.54;
%
% Text margin, uicontrol and axes gap (pixels):
mrg = 5;
gap = 4;
sid = 20;
%
% Slider position:
yid = 0;
ymx = 0;
%
% Handle of outlined text:
prv = [];
%
pmt = 'Click on a colorname...';
%
% Get palette names and colorspace functions:
[pnc,csf] = colornames();
%
if nargin<1
	idp = 1+rem(round(now*1e7),numel(pnc));
else
	palette = cnv1s2c(palette);
	assert(isChRo(palette),...
		'SC:colornames_view:palette:NotText',...
		'The first input <palette> must be a string scalar or a char row vector.')
	idp = strcmpi(palette,pnc);
	assert(any(idp),...
		'SC:colornames_view:palette:UnknownPalette',...
		'Palette "%s" is not supported. Call COLORNAMES() to list all palettes.',palette)
	idp = find(idp);
end
%
%% Color Sorting List %%
%
% For sorting:
ncs = {'AlphaNum';'Alphabet'};
ucs = {'Lab';'XYZ';'LCh';'YUV';'HSV';'RGB'};
% Get every permutation of the colorspaces:
cso = cellfun(@(s)perms(s(end:-1:1)),ucs,'UniformOutput',false);
cso = [ncs;cellstr(vertcat(cso{:}))];
acs = [ncs;ucs];
[lst,idl] = cellfun(@sort,acs,'UniformOutput',false);
%
if nargin<2
	ido = 1;
else
	order = cnv1s2c(order);
	assert(isChRo(order),...
		'SC:colornames_view:order:NotText',...
		'The second input <order> must be a string scalar or a char row vector.')
	ido = strcmpi(order,cso);
	assert(any(ido),...
		'SC:colornames_view:order:UnknownOption',...
		'The second input <order> must be one of:%s\b.',sprintf(' %s,',cso{:}))
	ido = find(ido);
end
%
% Intial color sorting index:
idx = 1:numel(colornames(pnc{idp}));
%
%% Create a New Figure %%
%
if isempty(fgh)||~ishghandle(fgh)
	% Figure with zoom and pan functions:
	fgh = figure('HandleVisibility','callback', 'IntegerHandle','off',...
		'NumberTitle','off', 'Name','ColorNames View', 'Color','white',...
		'Toolbar','figure', 'Units','pixels', 'Tag',mfilename, 'Visible','on');
	%
	fgp = get(fgh,'Position');
	inh = uicontrol(fgh, 'Units','pixels', 'Style','text', 'HitTest','off',...
		'Visible','on',	'String','Initializing the figure... please wait.');
	inx = get(inh,'Extent');
	set(inh,'Position',[fgp(3:4)/2-inx(3:4)/2,inx(3:4)])
	%
	% Axes and scrolling slider:
	slh = uicontrol('Parent',fgh, 'Style','slider', 'Visible','off',...
		'Enable','on', 'Value',1, 'Min',0, 'Max',1,...
		'FontUnits','pixels', 'Units','pixels', 'Callback',@cnvSldClBk);
	axh = axes('Parent',fgh, 'Visible','off', 'Units','pixels',...
		'YDir','reverse', 'XTick',[], 'YTick',[], 'XLim',[0,1], 'YLim',[0,1]);
	% Palette and color sorting method drop-down menus:
	pah = uicontrol('Parent',fgh, 'Style','popupmenu', 'String',pnc,...
		'ToolTip','Color Scheme', 'Units','pixels',...
		'Visible','off', 'Callback',@cnvPalClBk);
	srh = uicontrol('Parent',fgh, 'Style','popupmenu', 'String',cso,...
		'ToolTip','Sort Colors', 'Units','pixels',...
		'Visible','off', 'Callback',@cnvSrtClBk);
	coh = uicontrol('Parent',fgh, 'Style','edit', 'String',pmt,...
		'ToolTip','RGB Value',   'Units', 'pixels', 'Visible','off',...
		'HorizontalAlignment','left', 'Enable','inactive');
else
	set(coh,'String',pmt)
end
set(pah,'Value',idp);
set(srh,'Value',ido);
%
fgo = get(fgh, 'Pointer');
set(fgh, 'Pointer','watch')
drawnow()
%
%% Callback Functions %%
%
	function cnvPalClBk(h,~) % Palette Menu CallBack
		% Select a new palette.
		idp = get(h,'Value');
		set(slh, 'Value',1)
		set(coh, 'String',pmt)
		set(fgh, 'Pointer','watch')
		drawnow()    %disp('pal: drawing new text...')
		cnvTxtDraw() %disp('pal: sorting text...')
		cnvSortBy()  %disp('pal: resizing text...')
		cnvResize()  %disp('pal: complete!')
		set(fgh, 'Pointer',fgo)
	end
%
	function cnvSrtClBk(h,~) % Sort-Order Menu CallBack
		% Select the color sorting method.
		ido = get(h,'Value');
		set(fgh, 'Pointer','watch')
		drawnow()
		cnvSortBy()
		cnvResize()
		set(fgh, 'Pointer',fgo)
	end
%
	function cnvSldClBk(h,~) % Slider CallBack
		% Scroll the axes by changing the axes limits.
		tmp = diff(get(axh,'Ylim'));
		set(axh, 'Ylim',[0,tmp] + (ymx-tmp)*(1-get(h,'Value')));
	end
%
	function cnvZoomClBk(~,~) % Zoom CallBack
		% Change the font and margin sizes.
		tmp = diff(get(axh,'Ylim'));
		set(txh, 'FontSize',txs/tmp);
		set(txh, 'Margin',mrg/tmp);
	end
%
	function cnvPanClBk(~,~) % Pan CallBack
		% Move the scroll-bar to match panning of the axes.
		tmp = get(axh,'Ylim');
		set(slh, 'Value',max(0,min(1,1-tmp(1)/(ymx-diff(tmp)))))
	end
%
%% Color Sorting %%
%
	function cnvSortBy()
		[tmp,ids] = sort(cso{ido});
		idc = strcmp(tmp,lst);
		ids = ids(idl{idc});
		switch acs{idc}
			case 'AlphaNum'
				idx = csf.natsort(cnc);
				return
			case 'Alphabet'
				rxi = '^([-+]?\d+\w*)\s+(.+)$';
				tkn = regexp(cnc(:),rxi,'once','tokens');
				if all(cellfun('length',tkn)==2)
					tkn = vertcat(tkn{:});
					[~,idx] = sort(lower(tkn(:,2)));
				else
					[~,idx] = sort(lower(cnc(:)));
				end
				return
			case 'RGB'
				mat = rgb;
			case 'HSV'
				mat = csf.rgb2hsv(rgb);
			case 'XYZ'
				mat = csf.rgb2xyz(rgb);
			case 'Lab'
				mat = csf.xyz2lab(csf.rgb2xyz(rgb));
			case 'LCh'
				mat = csf.lab2lch(csf.xyz2lab(csf.rgb2xyz(rgb)));
			case 'YUV' % BT.709
				mat = csf.invgamma(rgb) * [...
					+0.2126, -0.19991, +0.61500;...
					+0.7152, -0.33609, -0.55861;...
					+0.0722, +0.43600, -0.05639];
			otherwise
				error('SC:colornames_view:space:UnknownOption',...
					'Colorspace "%s" is not supported.',cso{ido})
		end
		[~,idx] = sortrows(mat,ids);
	end
%
%% Re/Draw Text Strings %%
%
	function txtClBk(hnd,~,s,bgd)
		uistack(hnd,'top')
		try %#ok<TRYNC>
			set(prv, 'EdgeColor','none')
		end % setting one object is much faster than setting all objects.
		prv = hnd;
		%
		[R,G,B] = ndgrid(0:1);
		cgd = [R(:),G(:),B(:)];
		[~,idc] = max(sum(bsxfun(@minus,cgd,bgd).^2,2));
		hxs = sprintf('%02X',round(bgd*255));
		dcs = sprintf(',%.5f',bgd);
		set(coh, 'String',sprintf('#%s [%s] %s',hxs,dcs(2:end),s));
		set(hnd, 'EdgeColor',cgd(idc,:))
	end
%
txf = @(s,b,c)text('Parent',axh, 'String',s, 'BackgroundColor',b,...
	'Color',c,'Margin',mrg, 'Units','data', 'Interpreter','none',...
	'VerticalAlignment','bottom', 'HorizontalAlignment','right',...
	'Clipping','on', 'ButtonDownFcn',{@txtClBk,s,b},'LineWidth',3);
%
	function cnvTxtDraw()
		% Delete any existing colors:
		delete(txh(ishghandle(txh)))
		drawnow()
		% Get new colors:
		[cnc,rgb] = colornames(pnc{idp});
		% Calculate the text color:
		baw = (rgb*[0.298936;0.587043;0.114021])<thr;
		% Draw new colors in the axes:
		txh = cellfun(txf,cnc,num2cell(rgb,2),num2cell(baw(:,[1,1,1]),2),'Uni',0);
		txh = reshape([txh{:}],[],1);
		txs = get(txh(1),'FontSize');
		prv = txh(1);
	end
%
%% Resize the Axes and UIControls, Move the Colors %%
%
	function cnvResize(~,~)
		%
		%disp('rsz: preprocessing...')
		%
		zoom(fgh,'out');
		%
		if nargin
			set(fgh, 'Pointer','watch')
			drawnow()
		end
		%
		ecv = get(prv, 'EdgeColor');
		set(prv, 'EdgeColor','none');
		%
		set(axh, 'Ylim',[0,1])
		set(txh, 'Units','pixels', 'FontSize',txs, 'Margin',mrg)
		%disp('rsz: getting extents...')
		txe = cell2mat(get(txh(:),'Extent'));
		%disp('rsz: postprocessing...')
		top = get(slh,'FontSize')*2;
		fgp = get(fgh,'Position'); % [left bottom width height]
		hgt = round(fgp(4)-3*gap-top);
		wid = fgp(3)-3*gap-sid;
		pos = [gap,gap,wid,hgt];
		%
		% Calculate color lengths from text and margins:
		%disp('rsz: calculating text extent...')
		txw = 2*mrg+txe(idx,3);
		txc = cumsum(txw);
		%
		% Preallocate position array:
		txm = mean(txw);
		out = zeros(ceil(1.1*[txc(end)/pos(3),pos(3)/txm]));
		% Split colors into lines that fit the axes width:
		idb = 1;
		idr = 0;
		tmp = 0;
		while idb<=numel(txc)
			idr = idr+1;
			idq = max([idb,find((txc-tmp)<=pos(3),1,'last')]);
			out(idr,1:1+idq-idb) = txc(idb:idq)-tmp;
			tmp = txc(idq);
			idb = idq+1;
		end
		%
		% Calculate X and Y positions for each color:
		%disp('rsz: calulating text X & Y positions...')
		[~,txy,txx] = find(out.');
		txy = txy(:);
		txx = txx(:);
		yid = txy(end);
		txy = txy*(2*mrg+max(txe(idx,4)));
		ymx = txy(end)/pos(4);
		%
		% Resize the scrollbar, adjust scroll steps:
		nwp = [2*gap+wid,gap,sid,hgt];
		if ymx>1
			set(slh, 'Position',nwp, 'Enable','on', 'Value',1,...
				'SliderStep',max(0,min(1,[0.5,2]/(yid*(ymx-1)/ymx))))
		else
			set(slh, 'Position',nwp, 'Enable','off')
		end
		%
		% Resize the axes and drop-down menus:
		%disp('rsz: resizing the axes...')
		set(axh, 'Position',pos)
		uiw = (fgp(3)-gap)/4-gap;
		txw = 2*uiw+gap;
		lhs = gap+(0:2)*(uiw+gap);
		bot = fgp(4)-top-gap;
		set(pah, 'Position',[lhs(1),bot,uiw,top])
		set(srh, 'Position',[lhs(2),bot,uiw,top])
		set(coh, 'Position',[lhs(3),bot,txw,top])
		% Move text strings to the correct positions:
		%disp('rsz: moving text...')
		arrayfun(@(h,x,y)set(h,'Position',[x,y]),txh(idx),txx-mrg,pos(4)-txy+mrg);
		set(txh, 'Units','data')
		%
		set(prv, 'EdgeColor',ecv);
		%
		if nargin
			set(fgh, 'Pointer',fgo)
		end
		drawnow()
		%
		%disp('rsz: complete!')
	end
%
%% Initialize the Figure %%
%
%              disp('ini: drawing new text...')
cnvTxtDraw()  %disp('ini: sorting new text...')
cnvSortBy()   %disp('ini: resizing new text...')
cnvResize()   %disp('ini: complete!')
set([pah,srh,coh,slh], 'Visible','on')
set(fgh, 'Pointer',fgo, 'ResizeFcn',@cnvResize)
set(zoom(fgh), 'ActionPostCallback',@cnvZoomClBk);
set(pan(fgh),  'ActionPostCallback',@cnvPanClBk);
try %#ok<TRYNC>
	delete(inh)
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%colornames_view
function arr = cnv1s2c(arr)
% If scalar string then extract the character vector, otherwise data is unchanged.
if isa(arr,'string') && isscalar(arr)
	arr = arr{1};
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cnv1s2c