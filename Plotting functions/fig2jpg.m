function fig2jpg
% Convert a file that was previously a fig to a jpg

% Prompt folder
start_direct = 'D:\shear layer';
direct = uigetdir(start_direct, 'Choose Source Image Directory');

% Collect files and number of files
fig_files = dir([direct '\*.fig']);
num_files = length(fig_files);

% Open .fig and save as jpgs for space
for i = 1:num_files
    h = openfig([direct '\' fig_files(i).name]);
    h.Position = [300 300 800 600];
    saveas(h, [direct '\' fig_files(i).name(1:end-4)], 'jpg');
    close all
end

