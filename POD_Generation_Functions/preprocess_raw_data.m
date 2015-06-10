function [x, y, u, v, u_scale, l_scale] = preprocess_raw_data(x, y, u, v, l_scale, ...
    u_scale_gen, non_dim, xy_units, flip, image_range, direct)

% get number of images
images = size(u,3);

% Apply any flips to get the images in the correct orientation
for i = 1:images
    [x_temp, y_temp, u(:,:,i), v(:,:,i)] = image_rotation(x, y, u(:,:,i), v(:,:,i), flip);
end

x = x_temp;
y = y_temp;

% Crop images if requested
if ~isempty(image_range)
    x = x(image_range(1):image_range(2), image_range(3):image_range(4));
    y = y(image_range(1):image_range(2), image_range(3):image_range(4));
    u = u(image_range(1):image_range(2), image_range(3):image_range(4), :);
    v = v(image_range(1):image_range(2), image_range(3):image_range(4), :);
end

% Scale velocity by the inlet fast side streamwise velocity
if isa(u_scale_gen, 'function_handle')
    u_scale = u_scale_gen(u, direct);
else
    u_scale = u_scale_gen;
end

% if requested make values non-dimensionalized by u_scale l_scale
if non_dim
    u = u./u_scale;
    v = v./u_scale;

    % Change x & y from mm to meters
    x = x/(l_scale);
    y = y/(l_scale);
    
    l_scale = 1;
    u_scale = 1;
end

% if coordinates are in millimeters convert to meters
if strcmp(xy_units, 'mm')
    x = x/1000;
    y = y/1000;
end
end