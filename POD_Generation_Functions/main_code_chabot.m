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
% MAIN_CODE_CHABOT(..., IMAGE_RANGE) IMAGE_RANGE provides limits for PIV 
% images, in order to exclude portions of the flow that are not useful, may
% improve results, looking to to improve this further later
% 
% MAIN_CODE_CHABOT(..., DIRECT) if a DIRECT is provided will skip file
% selection GUI and perform analysis on files in the directory.
format long g
close all
clc

%%%%%%%%%%%% ALTER THIS TO MAKE PORTABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd('C:\Users\John-Desktop\Documents\MATLAB\thesis stuff\');

%% Set up function for given inputs, need to adjust to work better
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
        % Fourth parameter select dump2work status
        num_images      = varargin{1};
        overwrite       = varargin{2};
        save_pod        = varargin{3};
        dump2work       = varargin{4};
        save_figures    = 'none';
        image_range     = [];
        direct          = '';
    case 5
        % Fifth parameter select save file format
        num_images      = varargin{1};
        overwrite       = varargin{2};
        save_pod        = varargin{3};
        dump2work       = varargin{4};
        save_figures    = varargin{5};
        image_range     = [];
        direct          = '';
    case 6
        % Sixth paramater select PIV image range
        num_images      = varargin{1};
        overwrite       = varargin{2};
        save_pod        = varargin{3};
        dump2work       = varargin{4};
        save_figures    = varargin{5};
        image_range     = varargin{6};
        direct          = '';
    case 7
        % Seventh paramaeter select data directory
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

% Set up parrallel pool
if isempty(gcp)
    parpool;
end

gcp();

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

% find boundaries of velocity image
[~, ~, ~, ~, ~, ~, bnd_idx] = boundary_check_chabot(x, y, mean_u);

% Calculate volume elements of the mesh
vol_frac = voln_piv2(x, y, bnd_idx);

% Find fluxating velocity of flow
flux_u = zeros(size(u));
flux_v = zeros(size(v));

for i = 1:size(u,3);
    flux_u(:,:,i) = u(:,:,i) - mean_u;
    flux_v(:,:,i) = v(:,:,i) - mean_v;
end
clear u v

% Create a stacked data matrix for u and v velocities
flux_u      = reshape(flux_u, data_points, num_images);
flux_v      = reshape(flux_v, data_points, num_images);
mean_u      = reshape(mean_u, data_points, 1);
mean_v      = reshape(mean_v, data_points, 1);
vol_frac    = reshape(vol_frac, data_points, 1);

data.xg = x;
data.yg = y;

%% Perform Proper Orthogonal Decomposition
covariance = cal_covariance_mat(flux_u, flux_v, vol_frac);
[pod_u, pod_v, lambda2, eig_func, cutoff] =  calc_eig_modes2(covariance, flux_u, flux_v); 

pod_u = regroup(pod_u, dimensions);
pod_v = regroup(pod_v, dimensions);

% Flip images so they all "face" the same way
sign_flip = zeros(1, num_modes);
increment = 0.0000000000001;

% Figure out sign
for i = 1:num_modes
    sign_flip(i) = sign(mean(sign(pod_v(:,i)./(pod_v(:,i) + increment))));
end

% Apply flip
for i = 1:num_modes    
    pod_u(:,:,i) = pod_u(:,:,i)*sign_flip(i);
    pod_v(:,:,i) = pod_v(:,:,i)*sign_flip(i);
end

% Calculate vorticity
pod_u = reshape(pod_u, data_points, cutoff);
pod_v = reshape(pod_v, data_points, cutoff);

% Currently only calculating vorticity modes for num_modes may need to do
% for cutoff
vorticity = calc_vorticity(data, pod_u, pod_v, dimensions, bnd_idx);

%% Setup for Plotting and Plotting
if num_modes > 40
    num_plot = 40;
else
    num_plot = num_modes;
end

% Plot pod modes
Plotsvd2(data, pod_u(:,1:num_plot), dimensions, 'u', lambda2, bnd_idx, direct, save_figures);
Plotsvd2(data, pod_v(:,1:num_plot), dimensions, 'v', lambda2, bnd_idx, direct, save_figures);
Plotsvd2(data, vorticity(:,1:num_plot), dimensions, 'vorticity', lambda2, bnd_idx, direct, save_figures);
%% Save / Dump variables

% Save variables relavent to Galerkin to .mat files
if save_pod == true
    save([direct '\POD Data\POD.mat'], 'x', 'y', 'flux_u', 'flux_v', 'bnd_idx', ...
        'dimensions', 'eig_func', 'lambda2', 'mean_u', 'mean_v', 'pod_u', ...
        'pod_v', 'vol_frac', 'vorticity', 'cutoff');
end

% If requested place relvent galerkin variables in workspace
if dump2work == true
    putvar(x2, y2, bnd_idx, dimensions, eig_func, lambda2, mean_u, mean_v, ...
        pod_u, pod_v, vol_frac);
end

% return format
format short g
end