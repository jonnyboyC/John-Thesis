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
if num_modes > 0
    if custom
        full_path = [full_path filesep 'modes_' num2str(num_modes) '_custom'];
    else
        full_path = [full_path filesep 'modes_' num2str(num_modes)];
    end
    if ~isdir(full_path)
        % Exit if mode folder isn't found
        error('Galerkin coefficients for %d modes have not been produced', ...
            num_modes);
    end
end

% Get .mat wildcard
wildcard = get_wildcard(run_num, direct, data_folder, num_modes);

% Look in provided directory for .mat files
files = dir([full_path filesep wildcard]);

% If none are found prompt user to run previous level code
if size(files,1) == 0
    error(['Run required code to generate' data ' data for run ' num2str(run_num)]); 

% If only one is found use that as the file location
elseif size(files, 1) == 1
    file_loc = [full_path filesep files(1).name];

% If more than one is found, either select the provided mat_name from
% input or prompt user to select mat
elseif size(files, 1) > 1
   cd(full_path);
    fprintf(1, 'Please choose specific .mat file\n');
    file_name = uigetfile({'*.mat'}, 'Choose .mat file');
    file_loc = [full_path file_name];
end

end

