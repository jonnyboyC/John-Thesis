function [x, y, u, v] = Velocity_Read_Save(num_images, load_raw, image_range, flip, direct)
% VELOCITY_READ_PLOT_SAVE read raw num_images number of images from a selected
% directory in either .dat .vc7 or .mat formats
%
%   [x, y, u, v, num_x, num_y] = VELOCITY_READ_PLOT_SAVE(num_images)
% 

% Check to see if a saved file exists
if exist([direct filesep 'Processed Data' filesep 'Processed.mat'], 'file') == 2
    data = load([direct filesep 'Processed Data' filesep 'Processed.mat'], 'num_processed');
    num_processed = data.num_processed;
    if (load_raw == false && num_processed == num_images)
        load([direct filesep 'Processed Data' filesep 'Processed.mat'], 'x', 'y', 'u', 'v', 'num_x', 'num_y', 'u_scale');
        return;
    end
end

% Check the now set up folders for data
img_files = dir([direct filesep 'Raw Data' filesep '*']);

% Remove any directories from results
img_files = img_files([img_files.isdir]==0);
num_files = length(img_files);

% Get raw data file extension
[~, ~, file_type] = fileparts(img_files(1).name);

if num_images < num_files
    num_files = num_images;
end

% Load variables using appropriate method
if strcmpi(file_type, '.mat')
    [x, y, u, v] = load_mat(img_files, num_files, num_images, ...
                                    image_range, flip, direct);
elseif any(strcmpi(file_type, {'.vc7', '.VC7', '.im7', '.IM7'}))
    [x, y, u, v] = load_vc7(img_files, num_files, num_images, ...
                                    image_range, flip, direct);
else
    [x, y, u, v] = load_dat(img_files, num_files, num_images, ...
                                    image_range, flip, direct);
end
end

%% Load / Helper Functions

% Function to load files of .mat format
function [x, y, u, v] = load_mat(img_files, num_files, num_images, ...
                                        image_range, flip, direct)

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

% Function to load files of the .vc7/.im7 format
function [x, y, u, v] = load_vc7(img_files, num_files, num_images, ...
                                        image_range, flip, direct)

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

% Function to load files of the .dat in the format of cavity flow
function [x, y, u, v] = load_dat(img_files, num_files, num_images, ...
                                            image_range, flip, direct)

data_file = fopen([direct filesep 'Raw Data' filesep img_files(1).name]);

fgets(data_file); fgets(data_file);
fscanf(data_file,'%c',7);
num_x_org = fscanf(data_file,'%u',2);
num_x = num_x_org;

fscanf(data_file,'%c',4);
num_y_org = fscanf(data_file,'%u',2);
num_y = num_y_org;

fclose(data_file);


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
    data_file = fopen([direct filesep 'Raw Data' filesep file_name]);
    fgets(data_file); fgets(data_file); fgets(data_file);
    data = fscanf(data_file,'%g %g %g %g', [4 inf]);
    fclose(data_file);
    
    % Original file also takes the any additional files that contain a
    % * concatentated onto the end; such as B00001.vc7*. It also took
    % the file absolute path, currently these are not included

    xi = reshape(data(1,:), num_x_org, num_y_org);
    yi = reshape(data(2,:), num_x_org, num_y_org);
    ui = reshape(data(3,:), num_x_org, num_y_org);
    vi = reshape(data(4,:), num_x_org, num_y_org);  

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


% extract file name and number print current file
function [file_name, img_num] = update_progress(img_file)
    file_name = img_file.name;
    img_num = str2double(file_name(2:6));
    fprintf(1,'Processing file %s\n', file_name);
end