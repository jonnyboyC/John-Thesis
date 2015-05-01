function res = POD_Gen(varargin)
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
% problem.l_scale = 1
% Specify a scalar to be used as the length scale for the problem
%
% problem.u_scale_gen = 1
% Specify a scalar for function handle to calculate u_scale
%
% problem.flip = [false, false, false, false]
% Specify which directions and velocities are inverted
% [x, y, u, v]
%
% problem.update_bnds = false;
% If true launch gui modify boundaries
%
% problem.num_clusters = 10;
% set the number of clusters that should be determined
% 
% problem.exp_sampling_rate = 5;
% provide experimental sampling rate of the experiment, default is 5Hz

format long g
close all
clc

% List of fields that will be checked
fields = {  'num_images',   'load_raw',     'save_pod', ...
            'image_range',  'direct',       'l_scale', ...
            'u_scale_gen',  'save_figures', 'flip',...
            'update_bnds',  'num_clusters', 'exp_sampling_rate'};

% Parse problem structure provided to set it up correctly
if nargin == 1
    problem = parse_inputs(fields, @setdefaults_pod, varargin{1});
else
    fprintf('Provide a single structure as input, use help POD_Gen for information.\n');
    fprintf('Using Defaults\n\n');
    problem = parse_inputs(fields, @setdefaults_pod);
end

% Create more readable names
num_images  = problem.num_images;
load_raw    = problem.load_raw;
save_pod    = problem.save_pod;
image_range = problem.image_range;
direct      = problem.direct;
l_scale     = problem.l_scale;
u_scale_gen = problem.u_scale_gen;
save_figures= problem.save_figures;
flip        = problem.flip;
update_bnds = problem.update_bnds;
num_clusters= problem.num_clusters;
exp_sampling_rate = problem.exp_sampling_rate;

clear problem

% Check status of parrallel pool
if isempty(gcp('nocreate'));
    parpool('local');
end

%% Load and organize data

% Load simulation data from raw .vc7 or .mat, or from processed .mat
[x, y, u, v, u_scale, direct] = Velocity_Read_Save(num_images, load_raw, image_range, ...
                            l_scale, u_scale_gen, flip, direct);

% mean velocities and picture dimensions
mean_u = mean(u,3);
mean_v = mean(v,3);
dimensions = size(x);   

% Determine the number of images actually in memory
num_images = size(u,3);

% Determine resolution of velocity image
data_points = numel(x);

% Exactly define flow boundaries
[bnd_x, bnd_y, bnd_idx] = refine_bounds(x, y, mean_u, mean_v, direct, update_bnds);

% Calculate volume elements of the mesh
vol_frac = voln_piv2(x, y, bnd_idx);

% Check if mesh has even spacing
uniform = check_mesh(x, y);

% Find fluxating velocity of flow
flux_u = u - repmat(mean_u,1,1,num_images);
flux_v = v - repmat(mean_v,1,1,num_images);
    
clear u v

% Create a stacked data matrix for u and v velocities
flux_u      = reshape(flux_u, data_points, num_images);
flux_v      = reshape(flux_v, data_points, num_images);
mean_u      = reshape(mean_u, data_points, 1);
mean_v      = reshape(mean_v, data_points, 1);
vol_frac    = reshape(vol_frac, data_points, 1);

%% Perform Proper Orthogonal Decomposition

% Generate covariance matrix
covariance = cal_covariance_mat2(flux_u, flux_v, vol_frac, bnd_idx);

% Perform POD on fluctuating data
[pod_u, pod_v, lambda, modal_amp, cutoff] =  ...
    calc_eig_modes2(covariance, flux_u, flux_v); 

% Cluster resulting POD modes
[km_stoch, gm_stoch, gm_models, gm_groups, km_groups, centers] = ...
    cluster_POD(modal_amp, num_clusters, direct, save_figures);

pod_u = regroup(pod_u, dimensions);
pod_v = regroup(pod_v, dimensions);

% Flip images so they all "face" the same way
sign_flip = zeros(1, cutoff);

% Figure out sign and apply flip
for i = 1:cutoff
    sign_flip(i) = sign(mean(mean(sign(pod_v(:,:,i)./(pod_v(:,:,i) + eps)))));
    pod_u(:,:,i) = pod_u(:,:,i)*sign_flip(i);
    pod_v(:,:,i) = pod_v(:,:,i)*sign_flip(i);
end

% Calculate pod_vor
pod_u = reshape(pod_u, data_points, cutoff);
pod_v = reshape(pod_v, data_points, cutoff);

% Calculate voritcity
pod_vor = calc_pod_vor(pod_u, pod_v, dimensions, cutoff, bnd_idx, bnd_x, bnd_y, x, y);
pod_vor = reshape(pod_vor, data_points, cutoff);

%% Setup for Plotting and Plotting
if cutoff > 40
    num_plot = 40;
else
    num_plot = cutoff;
end

% Plot pod modes
if ~isempty(save_figures)
    Plotsvd2(x, y, pod_u(:,1:num_plot), dimensions, 'u', lambda, bnd_idx, direct, save_figures);
    Plotsvd2(x, y, pod_v(:,1:num_plot), dimensions, 'v', lambda, bnd_idx, direct, save_figures);
    Plotsvd2(x, y, pod_vor(:,1:num_plot), dimensions, 'vorticity', lambda, bnd_idx, direct, save_figures);
end

% Add mode zero
[modal_amp, lambda, pod_u, pod_v] = ...
    add_mode_zero(modal_amp, lambda, pod_u, pod_v, mean_u, mean_v);

%% Save / Return variables
run_num = floor(100000*rand(1));

% Place results in a structure
results.run_num = run_num;
results.exp_sampling_rate = exp_sampling_rate;

% mesh information
results.x = x;
results.y = y;
results.dimensions = dimensions;
results.bnd_idx = bnd_idx;
results.bnd_x = bnd_x;
results.bnd_y = bnd_y;
results.uniform = uniform;
results.vol_frac = vol_frac;

% fluctuating flow
results.flux_u = flux_u;
results.flux_v = flux_v;

% mean flow
results.mean_u = mean_u;
results.mean_v = mean_v;

% pod modes
results.pod_u = pod_u;
results.pod_v = pod_v;
results.pod_vor = pod_vor;

% k-mean clustering data
results.km_groups = km_groups;
results.centers = centers;
results.km_stoch = km_stoch;

% gaussian mixture model data
results.gm_models = gm_models;
results.gm_groups = gm_groups;
results.gm_stoch = gm_stoch;
results.cluster_range = 2:40;

% modal coordinates and energy captured, modes to 99%
results.modal_amp = modal_amp;
results.lambda = lambda;
results.cutoff = cutoff;

% scaling factoring
results.l_scale = l_scale;
results.u_scale = u_scale;


% Save variables relavent to Galerkin to .mat files
if save_pod == true
    save([direct filesep 'POD Data' filesep 'POD_run_' num2str(run_num) '.mat'], 'results', '-v7.3');
end

if nargout == 1
    res = results;
end

% return format
format short g
end