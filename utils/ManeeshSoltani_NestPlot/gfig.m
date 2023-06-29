function h = gfig(id, varargin)
% FIG - find or create a figure: [h=]fig('name')
%
% FIG and FIG(H) where H is a figure handle are almost identical to
% the built-in FIGURE, expect that if a Opt.new figure is created the
% name may include a Opt.prefix.
%
% FIG('name') looks for a figure with the specified name.  If it finds
% one it makes it current; if it fails it creates one with the given
% name and any Opt.prefix.
% 
% FIG(ID, ...) passes any unrecognised options to the figure.
%
% OPTIONS:
% 'Opt.new'		[off]	force Opt.new fig, appending sequence number to name
% 'Opt.prefix'	[hostname]	Opt.prefix for figure name
% See also: FIGURE.

% maneesh.
% pre-20010416: created
% 20130703: general Opt.prefix option; doc cleanup

% OPTIONS:
ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('new',[]);
ip.addParameter('prefix', 0);
ip.addParameter('figfunc', @sfigure);
ip.addParameter('grammObj', []);
ip.addParameter('clf', false);
ip.parse(varargin{:})
Opt = ip.Results;
figprops=util.struct2varargin(ip.Unmatched);

if isstring(id)
    id = char(id);
end

if nargin < 1
  h = figure;

  return;
end

if (~isempty(Opt.prefix) && Opt.prefix)
  if (~ischar(Opt.prefix))
    hostname=getenv('HOSTNAME');
    if(isempty(hostname))
      hostname=getenv('HOST');
    end
    if(isempty(hostname))
      [stat,hostname] = unix('uname -n');
    end
    hostname=strtok(hostname, '.');
    Opt.prefix = hostname;
  end
else
  Opt.prefix = 0;
end

if isstring(id)
    id = char(id);
end

if (~ischar(id))			% handle passed in
  if (ishandle(id))			% already a figure: leave it alone    
    h = Opt.figfunc(id);
  else					% set the Opt.new style name
    h = Opt.figfunc(id);
    if Opt.prefix
      set(h, 'name', sprintf('%s: %d', Opt.prefix, id));
    end
  end
else					% name passed in
  if Opt.prefix
    name = sprintf('%s: %s', Opt.prefix, id);
  else
    name = id;
  end

  h = findobj ('type', 'figure', 'name', name);
  if (~ isempty(h))			% already exists
    if (isempty(Opt.new))			% don't require a Opt.new figure
      set(0, 'CurrentFigure', h);	% doing it this way avoids focus issues
    else				% do require a Opt.new figure
      ii = 1;
      while 1				% look for a Opt.new name
        ii = ii + 1;
        name = sprintf('%s %d', id, ii);
	    h = findobj('type', 'figure', 'name', name);
        if isempty(h)
          h = Opt.figfunc('name', sprintf('%s: %s %d', hostname, id, ii));
          break;
        end
      end
    end
  else
    h = Opt.figfunc('name', name);
  end
end

if Opt.clf
    clf(h);
end
set(h, 'numbertitle', 'off', figprops{:});
titleFromName(h, Opt.grammObj);
