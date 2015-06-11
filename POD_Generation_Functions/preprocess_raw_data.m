function [X, U, u_scale, l_scale] = preprocess_raw_data(X, U, l_scale, ...
    u_scale_gen, non_dim, xy_units, flow_flip, image_range, direct)
% PREPROCESS_RAW_DATA perform the various scaling and image flipping, and
% cropping prior to actual calculation
%
%   [X, U, u_scale, l_scale] = preprocess_raw_data(X, U, l_scale,
%   u_scale_gen, non_dim, xy_units, flip, image_range, direct)

% Get fields
[x, u] = flow_comps(X, U);
comps = flow_ncomps(X);

% Apply any flips to get the images in the correct orientation
[X, U] = image_rotation(X, U, flow_flip);

% Crop images if requested
if ~isempty(image_range)
    % Quick check if crop index matches number of components
    if length(image_range) > comps
        error('provide a image_range vector of 2 * number of components');
    end
    
    % Get index range
    [X_idx, U_idx] = struct_index(image_range, 1:length(image_range), X, U);
    
    % Perform crops
    for i = 1:comps
        X.(x{i}) = X.(x{i})(X_idx{:});
        U.(u{i}) = U.(u{i})(U_idx{:});
    end
end

% Scale velocity by the inlet fast side streamwise velocity
if isa(u_scale_gen, 'function_handle')
    u_scale = u_scale_gen(U, direct);
else
    u_scale = u_scale_gen;
end

% if requested make values non-dimensionalized by u_scale l_scale
if non_dim
    for i = 1:comps
        X.(x{i}) = X.(x{i})./l_scale;
        U.(u{i}) = U.(u{i})./u_scale;
    end
    
    l_scale = 1;
    u_scale = 1;
end

% if coordinates are in millimeters convert to meters
if strcmp(xy_units, 'mm')
    for i = 1:comps
        X.(x{i}) = X.(x{i})./1000;
    end
end
end