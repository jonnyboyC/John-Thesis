function [pod_loc, direct] = prompt_folder()
    % Used by uigetdir to location initial folder    
    start_direct = 'D:\shear layer';
    
    % Prompt the user for location of Test folder
    fprintf(1, 'Please choose test data directory\n');
    direct = uigetdir(start_direct, 'Choose Source Image Directory'); 
    
    % Get file information 
    pod_files = dir([direct '\POD Data\*.mat']);
    
    % if no .mat are found prompt to run main_code_chabot
    if size(pod_files,1) == 0
       error('Run main_code_chabot to generated the POD data file'); 
       
    % if 1 .mat is found assume it is correct
    elseif size(pod_files, 1) == 1
        pod_loc = [direct '\POD Data\' pod_files(1).name];
    
    % IF 2 or more .mat are found prompt to check for correct .mat
    else
        error(['Multiple .mat files found in ' direct '\POD Data\ '...
            'check to ensure only one .mat exists'])
    end
end