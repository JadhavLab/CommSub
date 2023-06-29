function h = sfigure(varargin)
% SFIGURE  Create figure window (minus annoying focus-theft).
%
% Usage is identical to figure.
%
% Daniel Eaton, 2005
%
% See also figure

if nargin>=1 
	if ishandle(varargin{1})
		set(0, 'CurrentFigure', varargin{1});
    h=gcf;
	else
		h = figure(varargin{:});
	end
else
	h = figure;
end