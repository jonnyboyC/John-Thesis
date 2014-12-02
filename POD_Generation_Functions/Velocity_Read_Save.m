function [xi, yi, ui, vi, direct] = Velocity_Read_Save(num_images, overwrite, image_range, direct)
% VELOCITY_READ_PLOT_SAVE read num_images number of images from a selected
% directory.
%   [x, y, u, v, num_x, num_y] = VELOCITY_READ_PLOT_SAVE(num_images)
% 

%% Currently none of this used
% xls_file = 'A2_CpPHF_NS_20100728b.xls';
% xls_sheet = 'DataSheet';
% 
% % TODO have wind tunnel constants read from a file to make portable
% %[tun_constant tun_length, patm, num_phases, phase,count] = textread(...
% tunn_const  = 0.952;                % Tunnel constant
% len_scale   = 0.2032;               % Tunnel length meters
% tunn_pa     = 98600;                % Pa of wind tunnel
% phases      = 8;                    % Number of phases
% phase_cnt   = 0;                    % Counter for phases
% frame_dims  = [-70 -60 515 335];    % Frame size for all plotting
% 
% % read data from xls_file from sheet xls_sheet
% [data, ~]   = xlsread(xls_file, xls_sheet);
% 
% pres_cal    = 249;                  % Pressure calibration number
% acut_scale   = 1000;                 % scaling to get correct acutator position
% 
% run_num     = data(:,1);            % Chronological run number
% re          = data(:,2);            % Reynolds number by 1000 ie Re/1000
% angle_att   = data(:,3);            % Angle of attack
% dyn_pres    = pres_cal*data(:,4)/tunn_const;    % Dynamic Pressure
% stat_pres   = pres_cal*data(:,5);               % Static Pressure
% temp        = (data(:,6) + 459.67)/1.8;         % Temperature in K convert from F
% acut_x      = acut_scale*data(:,8); % Actuator Location
% pulse_vol   = data(:,11);           % pulser voltage
% pulse_freq  = data(:,12);           % pulser ferquency
% phase       = data(:,15);           % TODO phase ???
% 
% % throw away data we aren't using
% clear data
%% Currently in use

% Used by uigetdir to locate initial folder
start_direct = 'D:\shear layer';

% Prompt the user for location of Test Folder
fprintf(1, 'Please choose test data directory\n');
if strcmp(direct, '')
    direct = uigetdir(start_direct, 'Choose Source Image Directory');
end

% Set up the folders properly
update_folders(direct);

% change overwrtie of number of images requested changes
overwrite = update_overwrite(overwrite, num_images, direct);

% Check the now set up folders for data
saved_files = dir([direct '\Processed Data\*.mat']);
img_files = dir([direct '\Raw Data\*.vc7']);
num_files = length(img_files);
file_type = 'vc7';

% If files of the .vc7 are not found look for .im7 files
if num_files == 0
    img_files = dir([direct '\Raw Data\*.im7']);
    num_files = length(img_files);
    file_type = 'im7';
end 

% File no vc7 files are found look for .mat files
if num_files == 0
    img_files = dir([direct '\Raw Data\*.mat']);
    num_files = length(img_files);
    file_type = 'mat';
end

% Check to see if a saved file exists
if (overwrite == false && size(saved_files,1) == 2)
    load([direct '\Processed Data\' saved_files(2).name], 'xi', 'yi', 'ui', 'vi', 'num_x', 'num_y');
    return
end

if num_images < num_files
    num_files = num_images;
end

if strcmp(file_type, 'mat')
    [xi, yi, ui, vi] = load_mat(img_files, num_files, num_images, image_range, direct);
else
    [xi, yi, ui, vi] = load_vc7(img_files, num_files, num_images, image_range, direct);
end
end

% Function to load files of .mat format
function [xi, yi, ui, vi] = load_mat(img_files, num_files, num_images, image_range, direct)

% Get dimensions of image
load([direct '\Raw Data\' img_files(1).name], 'x', 'y');

if isempty(image_range)
    num_x = size(x,1);
    num_y = size(x,2);
else
    num_x = image_range(2)-image_range(1);
    num_y = image_range(4)-image_range(3);
end

% Preallocate Matrices
ui = zeros(num_x, num_y, num_files);
vi = zeros(num_x, num_y, num_files);

% Load images
for i = 1:num_files
   if ~img_files(i).isdir
       
        % Show current progress
        filename = update_progress(img_files(i));
        
        % Load individual images
        load([direct '\Raw Data\' filename], 'u', 'v');
        
        % Perform image rotation if necessary
        if isempty(image_range)
            [xi, yi, ui(:,:,i), vi(:,:,i)] = image_rotation(x, y, u, v);
        else
            [x_temp, y_temp, u_temp, v_temp] = image_rotation(x, y, u, v);
            xi = x_temp(image_range(1):image_range(2), image_range(3):image_range(4));
            yi = y_temp(image_range(1):image_range(2), image_range(3):image_range(4));
            ui(:,:,i) = u_temp(image_range(1):image_range(2), image_range(3):image_range(4),i);
            ui(:,:,i) = v_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        end
   end
end

% TODO LOOK into scaling of values

% Scale velocity by the inlet fast side streamwise velocity
u_scale = mean(ui,3);
u_scale = sort(u_scale(:));
u_scale = u_scale(floor(0.98*length(u_scale)));
ui = ui./u_scale;
vi = ui./v_scale;

% Change x & y from mm to meters
xi = xi/1000;
yi = yi/1000;

% Save Data to processed folder
num_processed = num_images;
save([direct '\Processed Data\Processed.mat'], 'xi', 'yi', 'ui', 'vi', 'num_x', 'num_y');
save([direct '\Processed Data\Num_Processed.mat'], 'num_processed');
end

% Function to load files of the .vc7/.im7 format
function [xi, yi, ui, vi] = load_vc7(img_files, num_files, num_images, image_range, direct)

% Get dimensions of image
lavdata = readimx([direct '\Raw Data\' img_files(1).name]);

if isempty(image_range)
    num_x = lavdata.Nx;
    num_y = lavdata.Ny;
else
    num_x = image_range(2)-image_range(1)+1;
    num_y = image_range(4)-image_range(3)+1;
end

% Preallocate matrices 
xi = zeros(num_x, num_y);  
yi = zeros(num_x, num_y);
ui = zeros(num_x, num_y, num_files);
vi = zeros(num_x, num_y, num_files);

% TODO original file was able to process multiple folders of data, will
% potentially want to add back in

% Load images
for i = 1:num_files
    % process if the file isn't a directory
    if ~img_files(i).isdir
        
        % Show current progress
        file_name = update_progress(img_files(i));
        
        % Original file also takes the any additional files that contain a
        % * concatentated onto the end; such as B00001.vc7*. It also took
        % the file absolute path, currently these are not included
        
        lavdata = readimx([direct '\Raw Data\' file_name]);
        [x,y,u,v] = showimx_mod(lavdata);

        % Rotate images to proper orientation
        if isempty(image_range)
            [xi, yi, ui(:,:,i), vi(:,:,i)] = image_rotation(x, y, u, v); 
        else
            [x_temp, y_temp, u_temp, v_temp] = image_rotation(x, y, u, v);
            xi = x_temp(image_range(1):image_range(2), image_range(3):image_range(4));
            yi = y_temp(image_range(1):image_range(2), image_range(3):image_range(4));
            ui(:,:,i) = u_temp(image_range(1):image_range(2), image_range(3):image_range(4));
            vi(:,:,i) = v_temp(image_range(1):image_range(2), image_range(3):image_range(4));
        end
    end
end

% TODO LOOK into scaling of values

% Scale velocity by the inlet fast side streamwise velocity
u_scale = mean(ui,3);
u_scale = sort(u_scale(:));
u_scale = u_scale(floor(0.98*length(u_scale)));
ui = ui./u_scale;
vi = vi./u_scale;

% Change x & y from mm to meters
xi = xi/1000;
yi = yi/1000;

% Save Data to processed folder
num_processed = num_images;
save([direct '\Processed Data\Processed.mat'], 'xi', 'yi', 'ui', 'vi', 'num_x', 'num_y');
save([direct '\Processed Data\Num_Processed.mat'], 'num_processed');

% TODO if we need to return to gathering data from multiple runs we should
% return u to a cell
end

% Check to see if number of images processed last time is same as this
% time. If it isn't update overwrite to true
function overwrite = update_overwrite(overwrite, num_images, direct)
    if overwrite == false
        if ~exist([direct '\Processed Data\Num_Processed.mat'], 'file')
            overwrite = true;
            return;
        end
        load([direct '\Processed Data\Num_Processed.mat'], 'num_processed');
        if num_images ~= num_processed
            overwrite = true;
        end
    end
end

% extract file name and number print current file
function [file_name, img_num] = update_progress(img_file)
    file_name = img_file.name;
    img_num = str2double(file_name(2:6));
    fprintf(1,'Processing file %s\n', file_name);
end