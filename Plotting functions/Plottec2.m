function [handle, cax] = Plottec2(data, h, bnd_idx)

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

if nargin == 1 || nargin == 3 || nargin == 4
    ax = newplot;
    h = surf(ax, data.x, data.y, zeros(size(data.x)), data.pod);
    h.FaceColor = 'interp';
    h.EdgeColor = 'none';
    if nargin == 3
        hold on
        scale = data.cmax;
        h2 = surf(ax, data.x, data.y, zeros(size(data.x)), scale*ones(size(data.pod)));
        hold off
        bnd_idx = double((bnd_idx == -1));
        h2.AlphaDataMapping = 'none';
        h2.FaceColor = 'flat';
        h2.EdgeColor = 'none';
        h2.AlphaData = bnd_idx;
        h2.FaceAlpha = 'flat';
    end
    set(ax, 'View', [0 90]);
    set(ax, 'Box', 'on');
    axis([minx maxx miny maxy]);
    colormap(ax, jet);
    axis(ax, 'equal')
    axis(ax, 'tight');
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