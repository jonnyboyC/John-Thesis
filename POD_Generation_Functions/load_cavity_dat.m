function [X, U] = load_cavity_dat(num_images, direct)
% LOAD_CAVITY_DAT load dat files in the format of the cavity into matlab 
%
%   [X, U] = LOAD_CAVITY_DAT(num_images, direct)
img_files = dir([direct filesep 'Raw Data' filesep '*']);

% Remove any directories from results
img_files = img_files([img_files.isdir]==0);
num_found = length(img_files);

if num_images < num_found
    num_found = num_images;
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

% Preallocate matrices 
x = zeros(num_x, num_y);  
y = zeros(num_x, num_y);
u = zeros(num_x, num_y, num_found);
v = zeros(num_x, num_y, num_found);

% Load images
for i = 1:num_found
    % Show current progress
    file_name = update_progress(img_files(i));
    
    % Open .dat file and read velocity snapshot
    data_file = fopen([direct filesep 'Raw Data' filesep file_name]);
    fgets(data_file); fgets(data_file); fgets(data_file);
    data = fscanf(data_file,'%g %g %g %g', [4 inf]);
    fclose(data_file);
    
    % Store each image into a matrix
    x = reshape(data(1,:), num_x_org, num_y_org);
    y = reshape(data(2,:), num_x_org, num_y_org);
    u(:,:,i) = reshape(data(3,:), num_x_org, num_y_org);
    v(:,:,i) = reshape(data(4,:), num_x_org, num_y_org);  
end   

% Place values into structure form
X.x = x;
X.y = y;

U.u = u;
U.v = v;
end