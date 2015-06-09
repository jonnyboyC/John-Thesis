function [x, y, u, v] = load_LaVision(img_files, num_files, num_images, ...
                                        image_range, flip, direct)
% LOAD_LAVISION load 2D .vc7 .im7 files into matlab, currently used for 
% mixing layer and airfoil data
%
%   [x, y, u, v] = LOAD_LAVISION(img_files, num_files, num_images,
%   image_range, flip, direct)  see help POD_GEN for information
                                    

% Get dimensions of image
lavdata = readimx([direct filesep 'Raw Data' filesep img_files(1).name]);

num_x = lavdata.Nx;
num_y = lavdata.Ny;

if ~isempty(image_range)
    if image_range(2) < num_x && image_range(1) >= 0
        num_x = image_range(2)-image_range(1)+1;
    end
    if image_range(4) < num_y && image_range(3) >= 0
        num_y = image_range(4)-image_range(3)+1;
    end
end

% Preallocate matrices 
x = zeros(num_x, num_y);  
y = zeros(num_x, num_y);
u = zeros(num_x, num_y, num_files);
v = zeros(num_x, num_y, num_files);

% Load images
for i = 1:num_files
    % Show current progress
    file_name = update_progress(img_files(i));

    % Original file also takes the any additional files that contain a
    % * concatentated onto the end; such as B00001.vc7*. It also took
    % the file absolute path, currently these are not included

    lavdata = readimx([direct filesep 'Raw Data' filesep file_name]);
    [xi,yi,ui,vi] = showimx_mod(lavdata);

    % Rotate images to proper orientation
    if isempty(image_range)
        [x, y, u(:,:,i), v(:,:,i)] = image_rotation(xi, yi, ui, vi, flip); 
    else
        [x_temp, y_temp, u_temp, v_temp] = image_rotation(xi, yi, ui, vi, flip);
        x = x_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        y = y_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        u(:,:,i) = u_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        v(:,:,i) = v_temp(image_range(1):image_range(2), image_range(3):image_range(4));
    end
end

% Save Data to processed folder
num_processed = num_images;
save([direct filesep 'Processed Data' filesep 'Processed.mat'], 'x', 'y', 'u', 'v', 'num_x', 'num_y', 'num_processed', '-v7.3');
end