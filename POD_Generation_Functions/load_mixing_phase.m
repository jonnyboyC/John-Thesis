function [X, U, num_processed] = load_mixing_phase(num_images, direct)
% LOAD_MIXING_PHASE load phase locked data from mat files into matlab, only
% used for data given by docter little, currently not used
%
%   [x, y, u, v] = LOAD_LAVISION(img_files, num_processed, num_images,
%   image_range, flip, direct) see help POD_GEN for information

% Check the now set up folders for data
img_files = dir([direct filesep 'Raw Data' filesep '*']);

% Remove any directories from results
img_files = img_files([img_files.isdir]==0);
num_processed = length(img_files);

if num_images < num_processed
    num_processed = num_images;
end

% Get dimensions of image
dims = load([direct filesep 'Raw Data' filesep img_files(1).name], 'x', 'y');
x = dims.x;
y = dims.y;

num_x = size(x,1);
num_y = size(x,2);

% Preallocate Matrices
u = zeros(num_x, num_y, num_processed);
v = zeros(num_x, num_y, num_processed);

% Load images
for i = 1:num_processed
    % Show current progress
    filename = update_progress(img_files(i));

    % Load individual images
    velocity = load([direct filesep 'Raw Data' filesep filename], 'u', 'v');
    u(:,:,i) = velocity.u;
    v(:,:,i) = velocity.v;
end

X.x = x;
X.y = y;

U.u = u;
U.v = v;
end