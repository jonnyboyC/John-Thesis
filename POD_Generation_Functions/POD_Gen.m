function results = POD_Gen(varargin)
% POD_GEN Creates and sets up proper folder format for vc7 images
% files, produces .mat files for processed vc7 files, and produced a
% requested number of POD modes with some simulaition specifications
%
% POD_GEN() prompts user for vc7 file path will process up to 1000
% images and producing figures of POD modes, without saving any information
% 
% POD_GEN(problem) Using fields provided in the stucture PROBLEM sets up
% simulation to specified by PROBLEM, all unfilled fields go to defaults
%
% problem has the following format with example values for default values
%
% problem.num_images = 1000 
% Specify the number of vc7 images or .mat matrices to be processed
%
% problem.load_raw = true 
% Should the previous processed mat be overwritten
%
% problem.save_pod = true
% Save the results of the simulation in POD.mat
%
% problem.save_figures = {'fig'}
% Request figures be saved as .fig or .jpg provided in a cell array
%
% problem.image_range = []
% Specify pixel range of images to be analyzed, in a matrix of form 
% [left, right, bottom, top]
%
% problem.direct = ''
% Specify absolute path of top of data directory. '' will prompt user for
% location
% 
% problem.l_scale = .3048
% Specify a scalar to be used as the lenght scale for the problem
%
% problem.u_scale = @u_scale_gen
% Specify a scalar for function handle to calculate u_scale

format long g
close all
clc

% List of fields that will be checked
fields = {  'num_images',   'load_raw',     'save_pod', ...
            'image_range',  'direct',       'l_scale', ...
            'u_scale_gen',  'save_figures' };

% Parse problem structure provided to set it up correctly
if nargin == 1
    problem = parse_inputs(fields, @setdefaults_pod, varargin{1});
else
    fprintf('Provide a single structure as input, use help POD_Gen for information.\n');
    fprintf('Using Defaults\n\n');
    problem = parse_inputs(fields, @setdefaults_pod);
end

% Create more convient names
num_images  = problem.num_images;
load_raw    = problem.load_raw;
save_pod    = problem.save_pod;
image_range = problem.image_range;
direct      = problem.direct;
l_scale     = problem.l_scale;
u_scale_gen = problem.u_scale_gen;
save_figures= problem.save_figures;

clear problem

% Set up parrallel pool
if isempty(gcp)
    parpool;
end

gcp();

%% Load and organize data

% Load simulation data from raw .vc7 or .mat, or from processed .mat
[x, y, u, v, u_scale, direct] = Velocity_Read_Save(num_images, load_raw, image_range, l_scale, u_scale_gen, direct);

% mean velocities and picture dimensions
mean_u = mean(u,3);
mean_v = mean(v,3);
dimensions = size(x);   

% Determine the number of images actually in memory
num_images = size(u,3);

% Determine resolution of velocity image
data_points = numel(x);

% find boundaries of velocity image
[~, ~, ~, ~, ~, ~, bnd_idx] = boundary_check_chabot(x, y, mean_u);

% Calculate volume elements of the mesh
vol_frac = voln_piv2(x, y, bnd_idx);

% Check if mesh has even spacing
uniform = check_mesh(x, y);

% Find fluxating velocity of flow
flux_u = zeros(size(u));
flux_v = zeros(size(v));

for i = 1:num_images;
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

data.x = x;
data.y = y;

%% Perform Proper Orthogonal Decomposition
covariance = cal_covariance_mat(flux_u, flux_v, vol_frac);
[pod_u, pod_v, lambda2, eig_func, cutoff] =  calc_eig_modes2(covariance, flux_u, flux_v); 

pod_u = regroup(pod_u, dimensions);
pod_v = regroup(pod_v, dimensions);

% Flip images so they all "face" the same way
sign_flip = zeros(1, cutoff);
increment = 0.0000000000001;

% Figure out sign and apply flip
for i = 1:cutoff
    sign_flip(i) = sign(mean(mean(sign(pod_v(:,:,i)./(pod_v(:,:,i) + increment)))));
    pod_u(:,:,i) = pod_u(:,:,i)*sign_flip(i);
    pod_v(:,:,i) = pod_v(:,:,i)*sign_flip(i);
end

vorticity = calc_vorticity2(pod_u, pod_v, dimensions, cutoff);

% Calculate vorticity
pod_u = reshape(pod_u, data_points, cutoff);
pod_v = reshape(pod_v, data_points, cutoff);
vorticity = reshape(vorticity, data_points, cutoff);

% Currently only calculating vorticity modes for num_modes may need to do
% for cutoff
% vorticity = calc_vorticity(data, pod_u, pod_v, dimensions, bnd_idx);

%% Setup for Plotting and Plotting
if cutoff > 40
    num_plot = 40;
else
    num_plot = cutoff;
end

% Plot pod modes
Plotsvd2(data, pod_u(:,1:num_plot), dimensions, 'u', lambda2, bnd_idx, direct, save_figures);
Plotsvd2(data, pod_v(:,1:num_plot), dimensions, 'v', lambda2, bnd_idx, direct, save_figures);
Plotsvd2(data, vorticity(:,1:num_plot), dimensions, 'vorticity', lambda2, bnd_idx, direct, save_figures);

%% Save / Dump variables
run_num = floor(100000*rand(1));

% Place results in a structure
results.run_num = run_num;
results.x = x;
results.y = y;
results.flux_u = flux_u;
results.flux_v = flux_v;
results.mean_u = mean_u;
results.mean_v = mean_v;
results.pod_u = pod_u;
results.pod_v = pod_v;
results.eig_func = eig_func;
results.lambda2 = lambda2;
results.l_scale = l_scale;
results.u_scale = u_scale;
results.vol_frac = vol_frac;
results.vorticty = vorticity;
results.cutoff = cutoff;
results.uniform = uniform;
results.dimensions = dimensions;
results.bnd_idx = bnd_idx;

% Save variables relavent to Galerkin to .mat files
if save_pod == true
    save([direct '\POD Data\POD.mat'], 'results', '-v7.3');
end

% return format
format short g
end