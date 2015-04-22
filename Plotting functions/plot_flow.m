function [h_surf, h_quiver, cax] = plot_flow(data, h_s, h_q)

% If no magnitude value is provided calculate the value
if ~isfield(data, 'pod') 
    if all(isfield(data, {'u', 'v'}))
        data.pod = sqrt(data.u.^2 + data.v.^2);
    else
        error('Must provide structure DATA with either u and v fields or pod field');
    end
end
if any(~isfield(data, {'x', 'y'}))
    error('Must provided structure DATA with both x and y fields');
end

% If a handle is provided simply update value instead of redrawling
if nargin == 1
    [h_s, ax] = Plottec2(data);
else
    h_s = Plottec2(data, h_s);
end

% Get no more than 50 quiver arrows
spacing_x = ceil(size(data.x, 1)/50);
spacing_y = ceil(size(data.x, 2)/50);

short_x = data.x(1:spacing_x:end, 1:spacing_y:end);
short_y = data.y(1:spacing_x:end, 1:spacing_y:end);
short_u = data.u(1:spacing_x:end, 1:spacing_y:end);
short_v = data.v(1:spacing_x:end, 1:spacing_y:end);

if nargin == 1 
    hold(ax, 'on')
    h_q = quiver(ax, short_x, short_y, short_u, short_v, 'color', [0 0 0]);
    hold(ax, 'off')
else
    h_q.UData = short_u;
    h_q.VData = short_v;
end

if isfield(data, 'format') && data.format == true
    ax.Title.String = 'Instanteous Flow Visualisation';
    ax.XLabel.String = 'x/D';
    ax.YLabel.String = 'y/D';
end

% Return quiver and surface handles
if nargout == 2
    h_surf = h_s;
    h_quiver = h_q;
end
if nargout == 3
    h_surf = h_s;
    h_quiver = h_q;
    cax = ax;
end
