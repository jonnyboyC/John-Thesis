function [X, U] = load_LaVision(num_images, direct)
% LOAD_LAVISION load 2D .vc7 .im7 files into matlab, currently used for 
% mixing layer and airfoil data
%
%   [x, y, u, v] = LOAD_LAVISION(img_files, num_processed, num_images,
%   image_range, flip, direct)  see help POD_GEN for information

% Check the now set up folders for data
img_files = dir([direct filesep 'Raw Data' filesep '*']);

% Remove any directories from results
img_files = img_files([img_files.isdir]==0);
num_found = length(img_files);

if num_images < num_found
    num_found = num_images;
end
                                  
% Get dimensions of image
lavdata = readimx([direct filesep 'Raw Data' filesep img_files(1).name]);

% Get mesh dimensions
num_x = lavdata.Nx;
num_y = lavdata.Ny;

% Preallocate matrices 
x = zeros(num_x, num_y);  
y = zeros(num_x, num_y);
u = zeros(num_x, num_y, num_found);
v = zeros(num_x, num_y, num_found);

% Load images
for i = 1:num_found
    
    % Show current progress
    file_name = update_progress(img_files(i));

    % Use proprietary readimx to get velocity data from im7 and bc7 files
    lavdata = readimx([direct filesep 'Raw Data' filesep file_name]);
    [x, y, u(:,:,i), v(:,:,i)] = showimx_mod(lavdata);
end

% Place quantities into structure
X.x = x;
X.y = y;

U.u = u;
U.v = v;
end