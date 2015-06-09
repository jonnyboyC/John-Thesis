function [x, y, u, v, num_processed] = load_cavity_dat(num_images, direct)
% LOAD_CAVITY_DAT load dat files in the format of the cavity into matlab 
%
%   [x, y, u, v] = LOAD_CAVITY_DAT(img_files, num_processed, num_images,
%   image_range, flip, direct) see help POD_GEN for information

% Check the now set up folders for data
img_files = dir([direct filesep 'Raw Data' filesep '*']);

% Remove any directories from results
img_files = img_files([img_files.isdir]==0);
num_processed = length(img_files);

if num_images < num_processed
    num_processed = num_images;
end

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
u = zeros(num_x, num_y, num_processed);
v = zeros(num_x, num_y, num_processed);

% Load images
for i = 1:num_processed
    % Show current progress
    file_name = update_progress(img_files(i));
    data_file = fopen([direct filesep 'Raw Data' filesep file_name]);
    fgets(data_file); fgets(data_file); fgets(data_file);
    data = fscanf(data_file,'%g %g %g %g', [4 inf]);
    fclose(data_file);
    
    % Original file also takes the any additional files that contain a
    % * concatentated onto the end; such as B00001.vc7*. It also took
    % the file absolute path, currently these are not included

    x = reshape(data(1,:), num_x_org, num_y_org);
    y = reshape(data(2,:), num_x_org, num_y_org);
    u(:,:,i) = reshape(data(3,:), num_x_org, num_y_org);
    v(:,:,i) = reshape(data(4,:), num_x_org, num_y_org);  
end   
end