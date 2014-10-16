function [xi, yi, ui, vi, direct] = Velocity_Read_Save(num_images, overwrite)
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
direct = uigetdir(start_direct, 'Choose Source Image Directory');

% Set up the folders properly
update_folders(direct);

% change overwrtie of number of images requested changes
overwrite = update_overwrite(overwrite, num_images, direct);

% Check the now set up folders for data
saved_files = dir([direct '\Processed Data\*.mat']);
img_files = dir([direct '\Raw Data\*.vc7']);
num_files = length(img_files);

% Check to see if a saved file exists
if (overwrite == false && size(saved_files,1) == 2)
    load([direct '\Processed Data\' saved_files(2).name], 'xi', 'yi', 'ui', 'vi', 'num_x', 'num_y');
    return
end

% If files of the .vc7 are not found look for .im7 files
if num_files == 0
    img_files = dir([direct '\Raw Data\*.im7']);
    num_files = length(img_files);
end 

if num_images < num_files
    num_files = num_images;
end

% Preallocat array size for number of image files
img_num = zeros(1,length(img_files));

lavdata = readimx([direct '\Raw Data\' img_files(1).name]);
num_x = lavdata.Nx;
num_y = lavdata.Ny;

xi = zeros(num_x, num_y, length(num_files));  
yi = zeros(num_x, num_y, length(num_files));
ui = zeros(num_x, num_y, length(num_files));
vi = zeros(num_x, num_y, length(num_files));

% TODO original file was able to process multiple folders of data, will
% potentially want to add back in

% TODO potentially to make parfor need to investigate more currently much
% faster than original impementation

for i = 1:num_files
    % process if the file isn't a directory
    if ~img_files(i).isdir
        
        % Show current progress
        [file_name, img_num(i)] = update_progress(img_files(i));
        
        % Original file also takes the any additional files that contain a
        % * concatentated onto the end; such as B00001.vc7*. It also took
        % the file absolute path, currently these are not included
        
        lavdata = readimx([direct '\Raw Data\' file_name]);
        [x,y,u,v] = showimx_mod(lavdata);
        % Need to check if this ever changes
        
        [xi(:,:,i), yi(:,:,i), ui(:,:,i), vi(:,:,i)] = image_rotation(x,y,u,v);            
    end
end

num_processed = num_images;
xi = xi(:,:,1);
yi = yi(:,:,1);
save([direct '\Processed Data\Processed.mat'], 'xi', 'yi', 'ui', 'vi', 'num_x', 'num_y');
save([direct '\Processed Data\Num_Processed.mat'], 'num_processed');

% TODO if we need to return to gathering data from multiple runs we should
% return u to a cell

end

% determine if the image is inverted in one or both of the axis
function [xn, yn, un, vn] = image_rotation(x, y, u, v)
    % prefill for speed
    xn = zeros(size(x));
    yn = zeros(size(y));
    un = zeros(size(u));
    vn = zeros(size(v));
    
    % create vector of indexes to avoid loops
    x_idx = 1:size(x,1);
    y_idx = 1:size(y,2);
    
    % perform flips
    if x(1,1) > x(end,1) && y(1,1) > y(1,end)
        xn(x_idx, y_idx) = x(end - x_idx + 1, end - y_idx + 1);
        yn(x_idx, y_idx) = y(end - x_idx + 1, end - y_idx + 1);
        un(x_idx, y_idx) = u(end - x_idx + 1, end - y_idx + 1);
        vn(x_idx, y_idx) = v(end - x_idx + 1, end - y_idx + 1);
    elseif x(1,1) > x(end,1) && y(1,1) < y(1, end)
        xn(x_idx,:) = x(end-x_idx+1,:);
        yn(x_idx,:) = y(end-x_idx+1,:);
        un(x_idx,:) = u(end-x_idx+1,:);
        vn(x_idx,:) = v(end-x_idx+1,:);
    elseif x(1,1) < x(end,1) && y(1,1) > y(1,end)
        xn(:,y_idx) = x(:,end-y_idx+1);
        yn(:,y_idx) = y(:,end-y_idx+1);
        un(:,y_idx) = u(:,end-y_idx+1);
        vn(:,y_idx) = v(:,end-y_idx+1);
    end
end

% Checks if the folders are properly set up 
function update_folders(direct)
    % Check if the data folders are present, create them if they are not.
    % These include a folder for raw vc7/similar files, a processed vc7
    % data, folder post pod calculation, and a figures folder
    if ~exist([direct '\Processed Data'], 'dir')
        mkdir(direct, '\Processed Data');
    end
    if ~exist([direct '\Raw Data'], 'dir')
        mkdir(direct, '\Raw Data');
    end
    if ~exist([direct '\POD Data'], 'dir')
        mkdir(direct, '\POD Data');
    end
    if ~exist([direct '\Figures'], 'dir')
       mkdir(direct, '\Figures');
       if ~exist([direct '\Figures\POD'], 'dir')
           mkdir(direct, '\Figures\POD');
       end
       if ~exist([direct '\Figures\Galerkin'], 'dir')
           mkdir(direct, '\Figures\Galerkin');
       end  
       if ~exist([direct '\Figures\Movies'], 'dir')
           mkdir(direct, '\Figures\Movies');
       end 
    end
    % Move any vc7 files into the Raw data folder
    [~, ~, ~] = movefile([direct '\*.vc7'], [direct '\Raw Data']); 
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