function [h_magnitude, h_direction, cax] = plot_vector_field(data, streamlines, h_mag, h_dir)
% PLOT_VECTOR_FIELD plot an individual vector field with quivers or
% streamlines to indicate direction
%
% [h_magnitude, h_direction, cax] = PLOT_VECTOR_FIELD(data, streamlines, h_mag, h_dir)
[x, u] = flow_comps_ip(data.X, data.U);
dims = flow_dims(data.X);

if dims ~= 2
    error('This function currently can only handle 2D meshes, need to revise for 3D mesh');
end
    
% If no magnitude value is provided calculate the value
if ~isfield(data, 'pod')
    
    %  Check that at least two components were passed in
    if isfield(data, {'U'}) && flow_ncomps(data.U) >= 2
        
        % Find the components that are in the plane of the image
        mag = 0;
        for i = 1:dims
            mag = mag + data.U.(u{i}).^2;
        end
        data.pod = sqrt(mag);
    else
        error('Must provide structure DATA with either u and v fields or pod field');
    end
end

% If a handle is provided simply update value instead of redrawling
if nargin == 2
    [h_mag, ax] = plot_scalar_field(data);
else
    h_mag = plot_scalar_field(data, h_mag);
end

% TODO STILL NEED TO CHECK ABOUT 3D MESH CASE HERE WE'RE JUST ASSUMING A 2D
% MESH

if streamlines
    if nargin == 2 
        h_dir = streamslice(data.X.(x{1})', data.X.(x{2})', data.U.(u{1})', data.U.(u{2})');
        set(h_dir,'Color','white')
    else
        delete(h_dir{1})
        h_dir = streamslice(data.X.(x{1})', data.X.(x{2})', data.U.(u{1})', data.U.(u{2})');
        set(h_dir,'Color','white')
    end
    h_dir = {h_dir};
else
    % Get no more than 50 quiver arrows
    spacing_x = ceil(size(data.X.(x{1}), 1)/50);
    spacing_y = ceil(size(data.X.(x{1}), 2)/50);
    
    short_idx = flow_index({[1 spacing_x 0], [1 spacing_y 0]}, [1, 2], data.X);
    
    for i = 1:dims
        short_X.(x{i}) = data.X.(x{i})(short_idx{:});
        short_U.(u{i}) = data.U.(u{i})(short_idx{:});
    end

    if nargin == 2 
        hold(ax, 'on')
        h_dir = quiver(ax, short_X.(x{1}), short_X.(x{2}),  short_U.(u{1}),  ...
            short_U.(u{2}), 'color', [1 1 1]);
        hold(ax, 'off')
    else
        h_dir.UData = short_U.(u{1});
        h_dir.VData = short_U.(u{2});
    end
end

if isfield(data, 'format') && data.format == true
    ax.Title.String = 'Instanteous Flow Visualisation';
    ax.XLabel.String = x(1);
    ax.YLabel.String = x(2);
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
