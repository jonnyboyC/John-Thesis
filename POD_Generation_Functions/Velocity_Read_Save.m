function [x, y, u, v] = Velocity_Read_Save(num_images, load_raw, image_range, flip, load_handle, direct)
% VELOCITY_READ_PLOT_SAVE read raw num_images number of images from a selected
% directory in either .dat .vc7 or .mat formats
%
%   [x, y, u, v] = VELOCITY_READ_PLOT_SAVE(num_images, load_raw, image_range, flip, direct)
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




