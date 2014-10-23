function [pod_loc, direct] = prompt_folder(data)
    data_folder = cell(length(data, 1));
    for i = 1:length(data)
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
    % Used by uigetdir to location initial folder    
    start_direct = 'D:\shear layer';
    
    % Prompt the user for location of Test folder
    fprintf(1, 'Please choose test data directory\n');
    direct = uigetdir(start_direct, 'Choose Source Image Directory'); 
    
    % Get file(s) information
    for i = 1:lenght(data);
        data_temp = check_data(data, i);
        files = dir([direct data_folder{i} '*.mat']);
    
        % if no .mat are found prompt to run main_code_chabot
        if size(files,1) == 0
            error(['Run required code to generate' data_temp ' data']); 
       
        % if 1 .mat is found assume it is correct
        elseif size(files, 1) == 1
            pod_loc = [direct data_folder files(1).name];
    
        % IF 2 or more .mat are found prompt to check for correct .mat
        else
            error(['Multiple .mat files found in ' direct data_folder{i}...
                ' check to ensure only one .mat exists'])
        end
    end
end

% Check to see if cell array
function data_temp = check_data(data, i)
if isCell(data)
    data_temp = data{i};
else
    data_temp = data;
end
end