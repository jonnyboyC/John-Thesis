function [file_name, img_num] = update_progress(img_file)
% UPDATE_PROGRES dispaly the current file number being read
%
%   [file_name, img_num] = UPDATE_PROGRESS(img_file)

% extract file name and number print current file
file_name = img_file.name;
img_num = str2double(file_name(2:6));
fprintf(1,'Processing file %s\n', file_name);
end