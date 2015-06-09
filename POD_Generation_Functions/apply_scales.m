function [xi, yi, ui, vi, u_scale, l_scale] = apply_scales(xi, yi, ui, vi, l_scale, ...
    u_scale_gen, non_dim, xy_units, direct)

% Scale velocity by the inlet fast side streamwise velocity
if isa(u_scale_gen, 'function_handle')
    u_scale = u_scale_gen(ui, direct);
else
    u_scale = u_scale_gen;
end

% if requested make values non-dimensionalized by u_scale l_scale
if non_dim
    ui = ui./u_scale;
    vi = vi./u_scale;

    % Change x & y from mm to meters
    xi = xi/(l_scale);
    yi = yi/(l_scale);
    
    l_scale = 1;
    u_scale = 1;
end

% if coordinates are in millimeters convert to meters
if strcmp(xy_units, 'mm')
    xi = xi/1000;
    yi = yi/1000;
end
end