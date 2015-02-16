function [xi, yi, ui, vi, u_scale, direct] = Velocity_Read_Save(num_images, load_raw, image_range, ...
                                            l_scale, u_scale_gen, flip_x, flip_y, direct)
% VELOCITY_READ_PLOT_SAVE read num_images number of images from a selected
% directory.
%   [x, y, u, v, num_x, num_y] = VELOCITY_READ_PLOT_SAVE(num_images)
% 

% Used by uigetdir to locate initial folder
start_direct = 'D:\shear layer';

% Prompt the user for location of Test Folder
fprintf(1, 'Please choose test data directory\n');
if strcmp(direct, '')
    direct = uigetdir(start_direct, 'Choose Source Image Directory');
end

% Set up the folders properly
update_folders(direct);

% change overwrtie of number of images requested changes
load_raw = update_overwrite(load_raw, num_images, direct);

% Check the now set up folders for data
img_files = dir([direct '\Raw Data\*']);

% Remove any directories from results
img_files = img_files([img_files.isdir]==0);
num_files = length(img_files);

% Get raw data file extension
[~, ~, file_type] = fileparts(img_files(1).name);

% Check to see if a saved file exists
if exist([direct '\Processed Data\Processed.mat'], 'file') == 2
    data = load([direct '\Processed Data\Processed.mat'], 'num_processed');
    num_processed = data.num_processed;
    if (load_raw == false && num_processed == num_images)
        load([direct '\Processed Data\Processed.mat'], 'xi', 'yi', 'ui', 'vi', 'num_x', 'num_y', 'u_scale');
        return
    end
end

if num_images < num_files
    num_files = num_images;
end

if strcmpi(file_type, '.mat')
    [xi, yi, ui, vi, u_scale] = load_mat(img_files, num_files, num_images, image_range, l_scale, u_scale_gen, flip_x, flip_y, direct);
elseif any(strcmpi(file_type, {'.vc7', '.im7'}))
    [xi, yi, ui, vi, u_scale] = load_vc7(img_files, num_files, num_images, image_range, l_scale, u_scale_gen, flip_x, flip_y, direct);
else
    [xi, yi, ui, vi, u_scale] = load_dat(img_files, num_files, num_images, image_range, l_scale, u_scale_gen, flip_x, flip_y, direct);
end
end

%% Load Functions

% Function to load files of .mat format
function [xi, yi, ui, vi, u_scale] = load_mat(img_files, num_files, num_images, image_range, ...
                                    l_scale, u_scale_gen, flip_x, flip_y, direct)

% Get dimensions of image
load([direct '\Raw Data\' img_files(1).name], 'x', 'y');

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
ui = zeros(num_x, num_y, num_files);
vi = zeros(num_x, num_y, num_files);

% Load images
for i = 1:num_files
    % Show current progress
    filename = update_progress(img_files(i));

    % Load individual images
    load([direct '\Raw Data\' filename], 'u', 'v');

    % Perform image rotation if necessary
    if isempty(image_range)
        [xi, yi, ui(:,:,i), vi(:,:,i)] = image_rotation(x, y, u, v, flip_x, flip_y);
    else
        [x_temp, y_temp, u_temp, v_temp] = image_rotation(x, y, u, v, flip_x, flip_y);
        xi = x_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        yi = y_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        ui(:,:,i) = u_temp(image_range(1):image_range(2), image_range(3):image_range(4),i);
        ui(:,:,i) = v_temp(image_range(1):image_range(2), image_range(3):image_range(4));
    end
end

% Apply any scaling that is present and translate x/y
[xi, yi, ui, vi, u_scale] = apply_scales(xi, yi, ui, vi, l_scale, u_scale_gen, direct);
[xi, yi] = correct_dims(xi, yi);

% Save Data to processed folder
num_processed = num_images;
save([direct '\Processed Data\Processed.mat'], 'xi', 'yi', 'ui', 'vi', 'num_x', 'num_y', 'u_scale', 'num_processed', '-v7.3');
end

% Function to load files of the .vc7/.im7 format
function [xi, yi, ui, vi, u_scale] = load_vc7(img_files, num_files, num_images, image_range, ...
                                    l_scale, u_scale_gen, flip_x, flip_y, direct)

% Get dimensions of image
lavdata = readimx([direct '\Raw Data\' img_files(1).name]);

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
xi = zeros(num_x, num_y);  
yi = zeros(num_x, num_y);
ui = zeros(num_x, num_y, num_files);
vi = zeros(num_x, num_y, num_files);

% TODO original file was able to process multiple folders of data, will
% potentially want to add back in

% Load images
for i = 1:num_files
    % Show current progress
    file_name = update_progress(img_files(i));

    % Original file also takes the any additional files that contain a
    % * concatentated onto the end; such as B00001.vc7*. It also took
    % the file absolute path, currently these are not included

    lavdata = readimx([direct '\Raw Data\' file_name]);
    [x,y,u,v] = showimx_mod(lavdata);

    % Rotate images to proper orientation
    if isempty(image_range)
        [xi, yi, ui(:,:,i), vi(:,:,i)] = image_rotation(x, y, u, v, flip_x, flip_y); 
    else
        [x_temp, y_temp, u_temp, v_temp] = image_rotation(x, y, u, v, flip_x, flip_y);
        xi = x_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        yi = y_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        ui(:,:,i) = u_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        vi(:,:,i) = v_temp(image_range(1):image_range(2), image_range(3):image_range(4));
    end
end

% Apply any scaling that is present and translate x/y
[xi, yi, ui, vi, u_scale] = apply_scales(xi, yi, ui, vi, l_scale, u_scale_gen, direct);
[xi, yi] = correct_dims(xi, yi);

% Save Data to processed folder
num_processed = num_images;
save([direct '\Processed Data\Processed.mat'], 'xi', 'yi', 'ui', 'vi', 'num_x', 'num_y', 'u_scale', 'num_processed', '-v7.3');
end

%TODO currently stub 
function [xi, yi, ui, vi, u_scale] = load_dat(img_files, num_files, num_images, image_range, ...
                                    l_scale, u_scale_gen, flip_x, flip_y, direct)

data_file = fopen([direct '\Raw Data\' img_files(1).name]);

fgets(data_file); fgets(data_file);
fscanf(data_file,'%c',7);
num_x_org = fscanf(data_file,'%u',2);
num_x = num_x_org;

fscanf(data_file,'%c',4);
num_y_org = fscanf(data_file,'%u',2);
num_y = num_y_org;

fclose(data_file);


if ~isempty(image_range)edi
    if image_range(2) < num_x && image_range(1) >= 0
        num_x = image_range(2)-image_range(1)+1;
    end
    if image_range(4) < num_y && image_range(3) >= 0
        num_y = image_range(4)-image_range(3)+1;
    end
end

% Preallocate matrices 
xi = zeros(num_x, num_y);  
yi = zeros(num_x, num_y);
ui = zeros(num_x, num_y, num_files);
vi = zeros(num_x, num_y, num_files);

% Load images
for i = 1:num_files
    % Show current progress
    file_name = update_progress(img_files(i));
    data_file = fopen([direct '\Raw Data\' file_name]);
    fgets(data_file); fgets(data_file); fgets(data_file);
    data = fscanf(data_file,'%g %g %g %g', [4 inf]);
    fclose(data_file);
    
    % Original file also takes the any additional files that contain a
    % * concatentated onto the end; such as B00001.vc7*. It also took
    % the file absolute path, currently these are not included

    x = reshape(data(1,:), num_x_org, num_y_org);
    y = reshape(data(2,:), num_x_org, num_y_org);
    u = reshape(data(3,:), num_x_org, num_y_org);
    v = reshape(data(4,:), num_x_org, num_y_org);

    

    % Rotate images to proper orientation
    if isempty(image_range)
        [xi, yi, ui(:,:,i), vi(:,:,i)] = image_rotation(x, y, u, v, flip_x, flip_y); 
    else
        [x_temp, y_temp, u_temp, v_temp] = image_rotation(x, y, u, v, flip_x, flip_y);
        xi = x_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        yi = y_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        ui(:,:,i) = u_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        vi(:,:,i) = v_temp(image_range(1):image_range(2), image_range(3):image_range(4));
    end
end   

% Apply any scaling that is present and translate x/y
[xi, yi, ui, vi, u_scale] = apply_scales(xi, yi, ui, vi, l_scale, u_scale_gen, direct);
[xi, yi] = correct_dims(xi, yi);

% Save Data to processed folder
num_processed = num_images;
save([direct '\Processed Data\Processed.mat'], 'xi', 'yi', 'ui', 'vi', 'num_x', 'num_y', 'u_scale', 'num_processed', '-v7.3');

end

%% Helper Functions

% Apply supplied scales
function [xi, yi, ui, vi, u_scale] = apply_scales(xi, yi, ui, vi, l_scale, u_scale_gen, direct)

% Scale velocity by the inlet fast side streamwise velocity
if isa(u_scale_gen, 'function_handle')
    u_scale = u_scale_gen(ui, direct);
else
    u_scale = u_scale_gen;
end

ui = ui./u_scale;
vi = vi./u_scale;

% Change x & y from mm to meters
xi = xi/(1000*l_scale);
yi = yi/(1000*l_scale);

end

function [xi, yi] = correct_dims(xi, yi)
xi = xi + abs(min(min(xi)));
yi = yi + abs(min(min(yi)));
end

% Check to see if number of images processed last time is same as this
% time. If it isn't update load_raw to true
function load_raw = update_overwrite(load_raw, num_images, direct)
    if load_raw == false
        if ~exist([direct '\Processed Data\Num_Processed.mat'], 'file')
            load_raw = true;
            return;
        end
        load([direct '\Processed Data\Num_Processed.mat'], 'num_processed');
        if num_images ~= num_processed
            load_raw = true;
        end
    end
end

% extract file name and number print current file
function [file_name, img_num] = update_progress(img_file)
    file_name = img_file.name;
    img_num = str2double(file_name(2:6));
    fprintf(1,'Processing file %s\n', file_name);
end