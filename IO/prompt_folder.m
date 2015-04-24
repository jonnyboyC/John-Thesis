function [file_loc, direct] = prompt_folder(data, run_num, direct, num_modes)
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

% If only one one input given prompt for directory
if nargin == 2    
    % Prompt the user for location of Test folder
    fprintf(1, ['Please choose test data ' data_folder '\n']);
    direct = uigetdir('', data_folder); 
end

% Get file(s) information
file_loc = get_data(data_folder, data, direct, run_num, num_modes);
end

function file_loc = get_data(data_folder, data, direct, run_num, num_modes)

% Get .mat wildcard
wildcard = get_wild(run_num, direct, data_folder, num_modes);

% If num_modes requested, look in Galerkin folder for mode data
full_path = [direct filesep data_folder];
if num_modes > 0
    full_path = [full_path filesep 'modes_' num2str(num_modes)];
    if ~isdir(full_path)
        % Exit if mode folder isn't found
        error('Galerkin coefficients for %d modes have not been produced', ...
            num_modes);
    end
end

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

% TODO clean this up
function wildcard = get_wild(run_num, direct, data_folder, num_modes)
if ischar(run_num)
    if strcmp(run_num, 'first')
        % Find newest file
        full_path = [direct filesep data_folder];
        
        % if num_modes specify move down the file tree
        if num_modes > 0
            full_path = [full_path filesep 'modes_' num2str(num_modes)];
            if ~isdir(full_path)
                % exit if folder not found
               error('Galerkin coefficients for %d modes have not been produced', ...
                   num_modes);
            end
        end
        
        % short files by order
        files = dir(full_path);
        [~, idx] = sort([files.datenum], 2, 'descend');
        files = files(idx);
        matches_idx = 0;
        
        % Return newest file not folder
        for i = 1:length(idx)
            if ~files(i).isdir 
                if ~matches_idx(1)
                    wildcard = files(i).name;
                    return;
                else
                if i == matches_idx(1)
                    wildcard = files(i).name;
                    return;
                end
                end
            end
        end
        
        % if not direct match return generic wild card
        wildcard = '*.mat';
    end
elseif isscalar(run_num)
    if run_num > 0
        wildcard = ['*' num2str(run_num) '*.mat'];
    else
        wildcard = '*.mat';
    end
else
    wildcard = '*.mat';
end
end