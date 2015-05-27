function [h_magnitude, h_direction, cax] = plot_vector_field(data, streamlines, h_mag, h_dir)
% PLOT_VECTOR_FIELD plot an individual vector field with quivers or
% streamlines to indicate direction
%
% [h_magnitude, h_direction, cax] = PLOT_VECTOR_FIELD(data, streamlines, h_mag, h_dir)

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
if nargin == 2
    [h_mag, ax] = plot_scalar_field(data);
else
    h_mag = plot_scalar_field(data, h_mag);
end

if streamlines
    if nargin == 2
        h_dir = streamslice(data.x', data.y', data.u', data.v');
        set(h_dir,'Color','black')
    else
        delete(h_dir{1})
        h_dir = streamslice(data.x', data.y', data.u', data.v');
        set(h_dir,'Color','black')
    end
    h_dir = {h_dir};
else
    % Get no more than 50 quiver arrows
    spacing_x = ceil(size(data.x, 1)/50);
    spacing_y = ceil(size(data.x, 2)/50);

    short_x = data.x(1:spacing_x:end, 1:spacing_y:end);
    short_y = data.y(1:spacing_x:end, 1:spacing_y:end);
    short_u = data.u(1:spacing_x:end, 1:spacing_y:end);
    short_v = data.v(1:spacing_x:end, 1:spacing_y:end);

    if nargin == 2 
        hold(ax, 'on')
        h_dir = quiver(ax, short_x, short_y, short_u, short_v, 'color', [0 0 0]);
        hold(ax, 'off')
    else
        h_dir.UData = short_u;
        h_dir.VData = short_v;
    end
end


if isfield(data, 'format') && data.format == true
    ax.Title.String = 'Instanteous Flow Visualisation';
    ax.XLabel.String = 'x/D';
    ax.YLabel.String = 'y/D';
end

% Return quiver and surface handles
if nargout == 2
    h_magnitude = h_mag;
    h_direction = h_dir;
end
if nargout == 3
    h_magnitude = h_mag;
    h_direction = h_dir;
    cax = ax;
end
