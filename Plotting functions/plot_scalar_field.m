function [handle, cax] = Plottec2(data, h)
% Create pseudo color plot, for flow visualization with included boundaries

if nargin == 1 || all(isfield(data, {'bnd_x', 'bnd_y'}))
    ax = newplot;
    h = surf(ax, data.x, data.y, zeros(size(data.x)), data.pod);
    h.FaceColor = 'interp';
    h.EdgeColor = 'none';
    if isfield(data, 'bnd_idx')
        % Get boundary index
        bnd_idx = data.bnd_idx;
        black = zeros([size(data.x), 3]);
        
        % Plot boundary index
        hold(ax, 'on');
        h2 = surf(ax, data.x, data.y, zeros(size(data.x)), black);
        hold(ax, 'off');
        
        % Make all non-boundary location transparent
        bnd_idx = double((bnd_idx == -1));
        h2.AlphaDataMapping = 'none';
        h2.FaceColor = 'flat';
        h2.EdgeColor = 'none';
        h2.AlphaData = bnd_idx;
        h2.FaceAlpha = 'flat';
    end
    if all(isfield(data, {'bnd_x', 'bnd_y'}))
        % Get Open Flow boundaries
        bnd_x = data.bnd_x;
        bnd_y = data.bnd_y;
        red = zeros(size(data.pod,1), size(data.pod,2), 3);
        red(:,:,1) = ones(size(data.pod));
        
        % Plot open flow boundaries
        hold(ax, 'on');
        h3 = surf(ax, data.x, data.y, zeros(size(data.x)), red);
        hold(ax, 'off');
        
        % Make all non-boundary location transparent
        flow_boundary = double(bnd_x ~= 0 | bnd_y ~= 0);
        h3.AlphaDataMapping = 'none';
        h3.FaceColor = 'flat';
        h3.EdgeColor = 'none';
        h3.AlphaData = flow_boundary;
        h3.FaceAlpha = 'flat';
    end
    
    % Determine boundaries of image
    minx = min(min(data.x));
    maxx = max(max(data.x));
    miny = min(min(data.y));
    maxy = max(max(data.y));
    
    % Modify viewing
    ax.View = [0 90];
    ax.Box = 'on';
    ax.NextPlot = 'replacechildren';
    colormap(ax, 'jet');
    axis(ax, [minx maxx miny maxy]);
    axis(ax, 'equal')
    axis(ax, 'tight');
    axis(ax, 'manual');
else
    % If plot already present update values
    h.CData = data.pod;
end
drawnow;

% Return handles
if nargout == 1
    handle = h;
end
if nargout == 2
    handle = h;
    cax = ax;
end