function [handle, cax] = Plottec2(data, h, varargin)
% Plot individual psuedo pcolor graphes

% Again may need to include zones
% Also may not the data structure at all
minx = min(min(data.x));
maxx = max(max(data.x));

miny = min(min(data.y));
maxy = max(max(data.y));

if nargin == 1 || nargin >= 3 
    ax = newplot;
    h = surf(ax, data.x, data.y, zeros(size(data.x)), data.pod);
    h.FaceColor = 'interp';
    h.EdgeColor = 'none';
    if size(varargin, 2) > 0
        bnd_idx = varargin{1};
        hold on
        h2 = surf(ax, data.x, data.y, zeros(size(data.x)), zeros(size(data.pod,1), size(data.pod,2), 3));
        hold off
        bnd_idx = double((bnd_idx == -1));
        h2.AlphaDataMapping = 'none';
        h2.FaceColor = 'flat';
        h2.EdgeColor = 'none';
        h2.AlphaData = bnd_idx;
        h2.FaceAlpha = 'flat';
    end
    if size(varargin, 2) == 3
        bnd_x = varargin{2};
        bnd_y = varargin{3};
        hold on 
        red = zeros(size(data.pod,1), size(data.pod,2), 3);
        red(:,:,1) = ones(size(data.pod));
        h3 = surf(ax, data.x, data.y, zeros(size(data.x)), red);
        hold off
        flow_boundary = double(bnd_x ~= 0 | bnd_y ~= 0);
        h3.AlphaDataMapping = 'none';
        h3.FaceColor = 'flat';
        h3.EdgeColor = 'none';
        h3.AlphaData = flow_boundary;
        h3.FaceAlpha = 'flat';
    end
    set(ax, 'View', [0 90]);
    set(ax, 'Box', 'on');
    axis([minx maxx miny maxy]);
    colormap(ax, jet);
    axis(ax, 'equal')
    axis(ax, 'tight');
else
    h.CData = data.pod;
end
drawnow;

if nargout == 1
    handle = h;
end
if nargout == 2
    handle = h;
    cax = ax;
end