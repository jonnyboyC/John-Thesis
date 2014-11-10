function [file_loc, direct] = prompt_folder(data, direct)
    data_folder = cell(size(data,2),1);
    file_loc = cell(size(data,1),1);
    
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
    
    % If only one one input given prompt for directory
    if nargin == 1
        % Used by uigetdir to location initial folder    
        start_direct = 'D:\shear layer';

        % Prompt the user for location of Test folder
        fprintf(1, 'Please choose test data directory\n');
        direct = uigetdir(start_direct, 'Choose Source Image Directory'); 
    end
    
    % Get file(s) information
    for i = 1:size(data,2);
        data_temp = check_data(data, i);
        files = dir([direct data_folder{i} '*.mat']);
    
        % if no .mat are found prompt to run main_code_chabot
        if size(files,1) == 0
            error(['Run required code to generate' data_temp ' data']); 
       
        % if 1 .mat is found assume it is correct
        elseif size(files, 1) == 1
            file_loc{i} = [direct data_folder{i} files(1).name];
    
        % IF 2 or more .mat are found prompt to check for correct .mat
        else
            error(['Multiple .mat files found in ' direct data_folder{i}...
                ' check to ensure only one .mat exists'])
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