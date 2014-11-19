function main_code_chabot(varargin)
% MAIN_CODE_CHABOT Creates and sets up proper folder format for vc7 images
% files, produces .mat files for processed vc7 files, and produced a
% requested number of POD modes
%
% MAIN_CODE_CHABOT() prompts user for vc7 file path will process up to 1000
% images and will not overwrite previous saves
%
% MAIN_CODE_CHABOT(NUM_IMAGES) loads NUM_IMAGES number of vc7 images, will
% by default will not overwrite previous .mat files
%
% MAIN_CODE_CHABOT(NUM_IMAGES, OVERWRITE) loads NUM_IMAGES number of vc7
% images will overwrite previous .mat files if OVERWRITE = TRUE regardless 
% of contents (default is FALSE)
%
% MAIN_CODE_CHABOT(NUM_IMAGES, OVERWRITE, SAVE_POD) loads NUM_IMAGES
% number of vc7 image will overwrite previous .mat files if OVERWRITE  =
% TRUE regardless of contents. If SAVE_POD = TRUE will save data used in
% galerkin projection to .mat file (default is TRUE)
%
% MAIN_CODE_CHABOT(NUM_IMAGES, OVERWRITE, SAVE_POD, DUMP2WORK)loads 
% NUM_IMAGES number of vc7 image will overwrite previous .mat files if 
% OVERWRITE = TRUE regardless of contents. If SAVE_POD = TRUE will save 
% data used in galerkin projection to .mat file (default is TRUE). If
% DUMP2WORK = TRUE relavent variables to Galerkin method will be placed
% into the workspace
%
% MAIN_CODE_CHABOT(..., SAVE_FIGURES) in addition to above if SAVE_FIGURES
% = 'fig' will save images as matlab .fig files if SAVE_FIGURES = 'jpg'
% will save images as jpgs.
%


% If SAVE_FIGURES = TRUE will save .fig of each pod mode generated.
% Set format and set to correct directory
format long g
close all
clc

%%%%%%%%%%%% ALTER THIS TO MAKE PORTABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd('C:\Users\John-Desktop\Documents\MATLAB\thesis stuff\');

%% Set up function for given inputs
switch nargin 
    case 0
        % Default: read 1000 images, don't overwrite, and don't save figures
        num_images      = 1000;
        overwrite       = false;
        save_pod        = true;
        dump2work       = false;
        save_figures    = 'none';
        image_range     = [];
        direct          = '';
    case 1
        % First parameter select number of images to load
        num_images      = varargin{1};
        overwrite       = false;
        save_pod        = true;
        dump2work       = false;
        save_figures    = 'none';
        image_range     = [];
        direct          = '';
    case 2
        % Second parameter select overwrite status
        num_images      = varargin{1};
        overwrite       = varargin{2};
        save_pod        = true;
        dump2work       = false;
        save_figures    = 'none';
        image_range     = [];
        direct          = '';
    case 3
        % Third parameter select save_figure status
        num_images      = varargin{1};
        overwrite       = varargin{2};
        save_pod        = varargin{3};
        dump2work       = false;
        save_figures    = 'none';
        image_range     = [];
        direct          = '';
    case 4
        num_images      = varargin{1};
        overwrite       = varargin{2};
        save_pod        = varargin{3};
        dump2work       = varargin{4};
        save_figures    = 'none';
        image_range     = [];
        direct          = '';
    case 5
        num_images      = varargin{1};
        overwrite       = varargin{2};
        save_pod        = varargin{3};
        dump2work       = varargin{4};
        save_figures    = varargin{5};
        image_range     = [];
        direct          = '';
    case 6
        num_images      = varargin{1};
        overwrite       = varargin{2};
        save_pod        = varargin{3};
        dump2work       = varargin{4};
        save_figures    = varargin{5};
        image_range     = varargin{6};
        direct          = '';
    case 7
        num_images      = varargin{1};
        overwrite       = varargin{2};
        save_pod        = varargin{3};
        dump2work       = varargin{4};
        save_figures    = varargin{5};
        image_range     = varargin{6};
        direct          = varargin{7};
    otherwise
        error('Too many input arguments');
end
pool = gcp;

% Some numbers i don't know what they mean
% mu0 = 94.5
% md = 1.296
% dp = md*mu0^2

%% Load and organize data
% Load velocity images from data, will load from raw files if processing
% has not been done, will load from .mat file otherwise. Select true to
% overwrite previous .mat files
[x, y, u, v, direct] = Velocity_Read_Save(num_images, overwrite, image_range, direct);

%%%%%%%%%%%%%%%%%% Temporary for testing will delete
% load('r1140i20aB002x10v80f1250p0_20100728b_n4.mat');
% load('r1140i20aB002x10v0f0p0_20100728b_n1.mat');
% [xcorrect,ycorrect, ~]=airfoil_rotation(20);
% x=x+xcorrect;
% y=y+ycorrect;
% num_images = 1000;
% u = u(:,1:num_images);
% v = v(:,1:num_images);
% u = reshape(u, Nx, Ny, num_images);
% v = reshape(v, Nx, Ny, num_images);
% x = reshape(x, Nx, Ny, 1);
% y = reshape(y, Nx, Ny, 1);
% for i = 1:size(u,3)
%     [xn, yn, un(:,:,i), vn(:,:,i)] = image_rotation(x, y, u(:,:,i), v(:,:,i));
% end
% x = xn;
% y = yn;
% u = un;
% v = vn;
% clear xn yn un vn
% direct = 'D:\shear layer\PIVData\Old Data';
%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_u = mean(u,3);
mean_v = mean(v,3);
dimensions = size(x);   

% Determine the number of images actually in memory
num_images = size(u,3);

% Determine resolution of velocity image
data_points = numel(x);

% TODO may want to request total number of modes
if num_images < 100
    num_modes = num_images;
else
    num_modes = 100;
end

% Again may need to do this in a cell fashion in case we need to read in
% multiple data

[~, ~, ~, ~, ~, ~, bnd_idx] = boundary_check_chabot(x, y, mean_u);
% [~, ~, ~, ~, ~, ~, bnd_idx] = boundary_check(x, y, mean_u);

vol_frac = voln_piv2(x, y, bnd_idx);

% Find fluxating velocity of flow
flux_u = zeros(size(u));
flux_v = zeros(size(v));

for i = 1:size(u,3);
    flux_u(:,:,i) = u(:,:,i) - mean_u;
    flux_v(:,:,i) = v(:,:,i) - mean_v;
end

% Create a stacked data matrix for u and v velocities
u           = reshape(u, data_points, num_images);
v           = reshape(v, data_points, num_images);
flux_u      = reshape(flux_u, data_points, num_images);
flux_v      = reshape(flux_v, data_points, num_images);
mean_u      = reshape(mean_u, data_points, 1);
mean_v      = reshape(mean_v, data_points, 1);
vol_frac    = reshape(vol_frac, data_points, 1);

data.xg = x;
data.yg = y;
% TODO calculate POD modes only for data points that not in the boundary,
% may potetially increase accuracy. Currently calculated modes including
% these boundaries

%% Perform Proper Orthogonal Decomposition
covariance = cal_covariance_mat(flux_u, flux_v, vol_frac);
[pod_u, pod_v, lambda2, eig_func] =  calc_eig_modes2(covariance, num_modes, flux_u, flux_v); 

pod_u1 = regroup(pod_u, dimensions);
pod_v1 = regroup(pod_v, dimensions);

% Flip images so they all "face" the same way
sign_flip = zeros(1, num_modes);
increment = 0.0000000000001;
for i = 1:num_modes
    sign_flip(i) = sign(mean(sign(pod_v(:,i)./(pod_v(:,i) + increment))));
end

for i = 1:num_modes    
    pod_u1(:,:,i) = pod_u1(:,:,i)*sign_flip(i);
    pod_v1(:,:,i) = pod_v1(:,:,i)*sign_flip(i);
end

pod_u1 = reshape(pod_u1, data_points, num_modes);
pod_v1 = reshape(pod_v1, data_points, num_modes);
vorticity = calc_vorticity(data, pod_u1, pod_v1, dimensions, bnd_idx);


%% Setup for Plotting and Plotting
if num_modes > 40
    num_plot = 40;
else
    num_plot = num_modes;
end


% Plot pod modes
Plotsvd2(data, pod_u1(:,1:num_plot), dimensions, 'u', lambda2, bnd_idx, direct, save_figures);
Plotsvd2(data, pod_v1(:,1:num_plot), dimensions, 'v', lambda2, bnd_idx, direct, save_figures);
Plotsvd2(data, vorticity(:,1:num_plot), dimensions, 'vorticity', lambda2, bnd_idx, direct, save_figures);
%% Save / Dump variables
% Save variables relavent to Galerkin to .mat files
if save_pod == true
    save([direct '\POD Data\POD.mat'], 'x', 'y', 'flux_u', 'flux_v', 'bnd_idx', ...
        'dimensions', 'eig_func', 'lambda2', 'mean_u', 'mean_v', 'pod_u1', ...
        'pod_v1', 'vol_frac', 'vorticity');
end
% If requested place relvent galerkin variables in workspace
x2 = x;
y2 = y;
if dump2work == true
    putvar(x2, y2, bnd_idx, dimensions, eig_func, lambda2, mean_u, mean_v, ...
        pod_u1, pod_v1, vol_frac);
end

% Need to look into if grid rotation is needed at all
% Need to ask what voln_piv does

format short g
end