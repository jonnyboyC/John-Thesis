function [file_loc, direct] = prompt_folder(data, direct, mat_name)
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

% Preallocate if cell
if iscell(data)
    data_folder = cell(size(data,2),1);
    file_loc = cell(size(data,2),1);
end

% Select correct subfolder
for i = 1:size(data,2)
    data_temp = check_data(data, i);
    switch data_temp;
        case 'POD'
            data_folder{i} = '\POD Data\';
        case 'Raw'
            data_folder{i} = '\Raw Data\';
        case 'Galerkin'
            data_folder{i} = '\Galerkin Coeff\';
    end
end

% if mat_name is empty fill with empty cells
if nargin < 3
    mat_name = cell(size(data,2),1);
end

% If only one one input given prompt for directory
if nargin == 1
    % Used by uigetdir to location initial folder    
    start_direct = 'D:\shear layer\PIVData';
    
    % Prompt the user for location of Test folder
    fprintf(1, ['Please choose test data ' data_folder{i}(2:end-1) '\n']);
    direct = uigetdir(start_direct, data_folder{i}(2:end-1)); 
end

% Get file(s) information
for i = 1:size(data,2);
    file_loc{i} = get_data(data_folder{i}, direct, mat_name);
end
    
end

function file_loc = get_data(data_folder, direct, mat_name)

% Look in provided directory for .mat files
files = dir([direct data_folder '*.mat']);

% If none are found prompt user to run previous level code
if size(files,1) == 0
    error(['Run required code to generate' data_temp ' data']); 

% If only one is found use that as the file location
elseif size(files, 1) == 1
    file_loc = [direct data_folder files(1).name];

% If more than one is found, either select the provided mat_name from
% input or prompt user to select mat
elseif size(files, 1) > 1
    if isempty(mat_name{1})
        cd([direct data_folder]);
        fprintf(1, 'Please choose .mat file for Galerkin\n');
        mat_name = uigetfile({'*.mat'}, 'Choose .mat file');
        file_loc = [direct data_folder mat_name];
    else
        file_loc = [direct data_folder mat_name];
    end
end

end

% Check to see if cell array
function data_temp = check_data(data, i)

if iscell(data)
    data_temp = data{i};
else
    data_temp = data;
end
end