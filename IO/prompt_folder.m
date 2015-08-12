function [file_loc, direct] = prompt_folder(data, run_num, direct, num_modes, custom)
% TODO really change this
% PROMPT_FOLDER Handles all the IO for selecting the correct data, can
% either select data manually by selection or with additional input
% argument select the data without prompt
%
% PROMPT_FOLDER(DATA) provide a string for DATA of 'POD', 'RAW' or 
% 'Galerkin' to manually select folder data will be pulled from, if 
% multiple .mat files are found prompt user to select specific folder
%
% PROMPT_FOLDER(DATA, DIRECT) provide a string for DATA of 'POD', 'RAW' or 
% 'Galerkin' with provided DIRECT to automatically collect .mat files. If
% multiple files are provided with manually prompt
% 
% PROMPT_FOLDER(DATA, DIRECT, MAT_NAME) same as prior but provides mat name
% explicitly 

% Select correct subfolder
switch data;
    case 'POD'
        data_folder = 'POD Data';
    case 'Raw'
        data_folder = 'Raw Data';
    case 'Galerkin'
        data_folder = 'Galerkin Coeff';
end

% if mat_name is empty fill with empty cells
if nargin < 4
    num_modes = 0;
end

if nargin < 5
    custom = false;
end

% If only one one input given prompt for directory
if nargin == 2    
    % Prompt the user for location of Test folder
    fprintf(1, ['Please choose test data ' data_folder '\n']);
    direct = uigetdir('', data_folder); 
end

% Get file(s) information
file_loc = get_data(data_folder, data, direct, run_num, num_modes, custom);
end

function file_loc = get_data(data_folder, data, direct, run_num, num_modes, custom)

% If num_modes requested, look in Galerkin folder for mode data
full_path = [direct filesep data_folder];

% If number of modes select, ie galerkin or Mod POD find mode folder
if num_modes ~= 0
    full_path = get_mode_folder(full_path, num_modes, custom);
end

% Get .mat wildcard
file = get_wildcard(run_num, full_path);

% Look in provided directory for .mat files
file_path = dir([full_path filesep file]);

% If none are found prompt user to run previous level code
if size(file_path,1) == 0
    error(['Run required code to generate' data ' data for run ' num2str(run_num)]); 

% If only one is found use that as the file location
elseif size(file_path, 1) == 1
    file_loc = [full_path filesep file_path(1).name];

% If more than one is found, either select the provided mat_name from
% input or prompt user to select mat
elseif size(file_path, 1) > 1
   ocd(full_path);
    fprintf(1, 'Please choose specific .mat file\n');
    file_name = uigetfile({'*.mat'}, 'Choose .mat file');
    file_loc = [full_path file_name];
end

end

