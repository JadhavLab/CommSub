function nestresize(ax, ppos)
% NESTRESIZE - resize a nestgroup: nestresize(AX, ['POS'])
%
% NESTRESIZE is primarily intended as a callback to allow automatic
% resizing of nested plots when the container is resized.  

% maneesh.
%
% 20130704: created
% 20140515: corrected recursive call arglist 


nestdata = getappdata(ax, 'NestData');
if isempty(nestdata)
  error('axes are neither nest container nor nested');
end

if nestdata.nestable == 0
  % call recursively on container
  if nargin > 1
    nestresize(nestdata.container, ppos);
  else
    nestresize(nestdata.container);
  end
  return;
end

if nargin > 1
  % position passed in -- call set, which will call us recursively
  set(ax, 'position', ppos);
  return;
end
  
% called with only container handle -- probably via callback

ppos = get(ax, 'position');

for cc = nestdata.children(:)'
  if ~isempty(cc)
    cdata = getappdata(cc, 'NestData');
    npos = cdata.position;
    cpos = [ppos(1:2) + npos(1:2).*ppos(3:4) npos(3:4).*ppos(3:4)];
    set(cc, 'position', cpos);
  end
end
% Ryan --- conditioning on cc's emptiness instead
% for cc = nestdata.children(:)'
%   if cc
%     cdata = getappdata(cc, 'NestData');
%     npos = cdata.position;
%     cpos = [ppos(1:2) + npos(1:2).*ppos(3:4) npos(3:4).*ppos(3:4)];
%     set(cc, 'position', cpos);
%   end
% end
    



