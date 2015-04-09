function [handle, cax] = plot_flow(data, h, format)

if ~isfield(data, 'pod') || all(isfield(data, {'u', 'v'}))
    data.pod = sqrt(u.^2 + v.^2);
else
    error('Must provide structure DATA with either u and v fields or pod field');
end
if any(~isfield(data, {'x', 'y'}))
    error('Must provided structure DATA with both x and y fields');
end


[h, ax] = Plottec2(data, h);


spacing_x = ceil(size(x, 1)/50);
spacing_y = ceil(size(x, 2)/50);

short_x = data.x(1:spacing_x:end, 1:spacing_y:end);
short_y = data.y(1:spacing_x:end, 1:spacing_y:end);
short_u = data.u(1:spacing_x:end, 1:spacing_y:end, image);
short_v = data.v(1:spacing_x:end, 1:spacing_y:end, image);

hold(ax, 'on')
quiver(ax, short_x, short_y, short_u, short_v, 'color', [0 0 0]);
hold(ax, 'off')

if nargin == 3 && format == true
    ax.Title.String = ['Instanteous Flow Visualisation for Snapshot ' num2str(image)];
    ax.XLabel.String = 'x/D';
    ax.YLabel.String = 'y/D';
end

if nargout == 1
    handle = h;
end
if nargout == 2
    handle = h;
    cax = ax;
end
