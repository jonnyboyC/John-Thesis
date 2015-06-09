function [x, y, u, v] = load_mixing_phase(img_files, num_files, num_images, ...
                                        image_range, flip, direct)
% LOAD_MIXING_PHASE load phase locked data from mat files into matlab, only
% used for data given by docter little, currently not used
%
%   [x, y, u, v] = LOAD_LAVISION(img_files, num_files, num_images,
%   image_range, flip, direct) see help POD_GEN for information

% Get dimensions of image
dims = load([direct filesep 'Raw Data' filesep img_files(1).name], 'x', 'y');
x = dims.x;
y = dims.y;

num_x = size(x,1);
num_y = size(x,2);

if ~isempty(image_range)
    if image_range(2) < num_x && image_range(1) >= 0
        num_x = image_range(2)-image_range(1)+1;
    end
    if image_range(4) < num_y && image_range(3) >= 0
        num_y = image_range(4)-image_range(3)+1;
    end
end

% Preallocate Matrices
u = zeros(num_x, num_y, num_files);
v = zeros(num_x, num_y, num_files);

% Load images
for i = 1:num_files
    % Show current progress
    filename = update_progress(img_files(i));

    % Load individual images
    velocity = load([direct filesep 'Raw Data' filesep filename], 'u', 'v');
    ui = velocity.u;
    vi = velocity.v;
    
    % Perform image rotation if necessary
    if isempty(image_range)
        [x, y, u(:,:,i), v(:,:,i)] = image_rotation(xi, yi, ui, vi, flip);
    else
        [x_temp, y_temp, u_temp, v_temp] = image_rotation(xi, yi, ui, vi, flip);
        x = x_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        y = y_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        u(:,:,i) = u_temp(image_range(1):image_range(2), image_range(3):image_range(4),i);
        v(:,:,i) = v_temp(image_range(1):image_range(2), image_range(3):image_range(4));
    end
end

% Save Data to processed folder
num_processed = num_images; %#ok<*NASGU>
save([direct filesep 'Processed Data' filesep 'Processed.mat'], 'x', 'y', 'u', 'v', 'num_x', 'num_y', 'num_processed', '-v7.3');
end