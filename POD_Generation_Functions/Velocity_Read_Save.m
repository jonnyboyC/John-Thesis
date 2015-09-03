function [X, U] = Velocity_Read_Save(num_images, load_raw, load_handle, direct)
% VELOCITY_READ_PLOT_SAVE read raw num_images number of images from a selected
% directory in either .dat .vc7 or .mat formats
%
%   [X, U] = VELOCITY_READ_PLOT_SAVE(num_images, load_raw, load_handle,
%   direct) checks if velocity data has already been stored in a .mat
%   otherwise uses load_handle to load velocity data

%#ok<*NASGU>

% Check to see if a saved file exists
if exist([direct filesep 'Processed Data' filesep 'Processed.mat'], 'file') == 2
    data = load([direct filesep 'Processed Data' filesep 'Processed.mat'], 'num_processed');
    num_processed = data.num_processed;
    if (load_raw == false && num_processed == num_images)
        load([direct filesep 'Processed Data' filesep 'Processed.mat'], 'X', 'U');
        return;
    end
end

% Load variables using appropriate method provided by the function handle
[X, U] = load_handle(num_images, direct);
num_processed = num_images; 

% Save Data to processed folder
save([direct filesep 'Processed Data' filesep 'Processed.mat'], 'X', 'U', 'num_processed', '-v7.3');
end




