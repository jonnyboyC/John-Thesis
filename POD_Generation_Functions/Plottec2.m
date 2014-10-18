function [handle, cax] = Plottec2(data, h)

% Again may need to include zones
% Also may not the data structure at all
data.x=data.xg;
data.y=data.yg;

minx = min(min(data.x));
maxx = max(max(data.x));

miny = min(min(data.y));
maxy = max(max(data.y));

% TODO check what rmghost does, currently no trace of this
% if blah ==0;
%    data.y = rmghost(data.y);
%    data.x = rmghost(data.x);
% end

if nargin == 1
    ax = newplot;
    h = surf(data.x, data.y, data.pod);
    set(ax, 'View', [0 90]);
    set(ax, 'Box', 'on');
    axis([minx maxx miny maxy]);
    colormap(jet);
    shading interp
    axis equal
else
    h.ZData = data.pod;
end
drawnow;
    
if nargout == 1
    handle = h;
end

if nargout == 2
    handle = h;
    cax = ax;
end