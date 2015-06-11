function [res_pod, res_clust] = POD_Gen(varargin)
% POD_GEN Creates and sets up proper folder format for vc7 images
% files, produces .mat files for processed vc7 files, and produced a
% requested number of POD modes with some simulaition specifications
%
%   POD_GEN() prompts user for vc7 file path will process up to 1000
%   images and producing figures of POD modes, without saving any 
%   information
% 
%   POD_GEN(problem) Using fields provided in the stucture PROBLEM sets up
%   simulation to specified by PROBLEM, all unfilled fields go to defaults
%
%   problem has the following format with example values for default values
%
%   problem.num_images = 1000 
%   Specify the number of vc7 images or .mat matrices to be processed
%
%   problem.load_raw = true 
%   Should the previous processed mat be overwritten
%
%   problem.save_pod = true
%   Save the results of the simulation in POD.mat
%
%   problem.save_figures = {'fig'}
%   Request figures be saved as .fig, .jpg or .png provided in a cell array
%
%   problem.image_range = []
%   Specify pixel range of images to be analyzed, in a matrix of form 
%   [left, right, bottom, top]
%
%   problem.direct = ''
%   Specify absolute path of top of data directory. '' will prompt user for
%   location
% 
%   problem.cluster = true
%   Specify if you want to cluster results or not
%
%   problem.average_mesh = false
%   Specify if you would like to average raw data over 2x2 range
%
%   problem.l_scale = 1
%   Specify a scalar to be used as the length scale for the problem
%
%   problem.u_scale_gen = 1
%   Specify a scalar for function handle to calculate u_scale
%
%   problem.flip = {}
%   Specify which directions and velocities are inverted as a cell array
%   such as {'u', 'w', 'y'}
%
%   problem.update_bnds = false;
%   If true launch gui to modify boundaries
%
%   problem.filter = false;
%   If true apply a guide filter to images, to attempt to minimize affect of
%   PIV artifacts
%
%   problem.num_clusters = 10;
%   set the number of clusters that should be determined
% 
%   problem.exp_sampling_rate = 5;
%   provide experimental sampling rate of the experiment, default is 5Hz
%
%   problem.streamlines = true
%   plot mode using streamline instead of quivers
%
%   problem.non_dim = false
%   non-dimensionalize the units by u/u_scale x/l_scale
%
%   problem.xy_units = 'mm'
%   specify the units the x and y coordinate are in currently accepts 'm'
%   for meters and 'mm' for millimeters

format long g
%close all
clc

% List of fields that will be checked
fields = {  'num_images',   'load_raw',     'save_pod', ...
            'image_range',  'direct',       'l_scale', ...
            'u_scale_gen',  'save_figures', 'flow_flip',...
            'update_bnds',  'num_clusters', 'exp_sampling_rate',...
            'cluster',      'average_mesh', 'filter', ...
            'streamlines',  'non_dim',      'xy_units', ...
            'load_handle',  'open_flow'};

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
load_handle = problem.load_handle;
save_pod    = problem.save_pod;
image_range = problem.image_range;
direct      = problem.direct;
l_scale     = problem.l_scale;
u_scale_gen = problem.u_scale_gen;
save_figures= problem.save_figures;
flow_flip   = problem.flow_flip;
update_bnds = problem.update_bnds;
num_clusters= problem.num_clusters;
cluster     = problem.cluster;
filter      = problem.filter;
streamlines = problem.streamlines;
average_mesh= problem.average_mesh;
non_dim     = problem.non_dim;
xy_units    = problem.xy_units;
open_flow   = problem.open_flow;
exp_sampling_rate = problem.exp_sampling_rate;

clear problem

% Check status of parrallel pool
if isempty(gcp('nocreate'));
    parpool('local');
end

%% Load and Preprocess Data

% Prompt the user for location of folder if not provided
fprintf(1, 'Please choose test data directory\n');
if strcmp(direct, '')
    direct = uigetdir('', 'Choose Source Image Directory');
end

% Update folder to correct format
update_folders(direct);

% Load simulation data from raw .vc7 or .mat, or from processed .mat
[X, U] = Velocity_Read_Save(num_images, load_raw, load_handle, direct);  

% Apply scaling to flow variables
[X, U, u_scale, l_scale] = preprocess_raw_data(X, U, l_scale, u_scale_gen, ...
                                            non_dim, xy_units, flow_flip, image_range, direct);

% Check if mesh has even spacing
uniform = check_mesh(X);

if uniform == false
    streamlines = false;
end

% If request compress mesh by a factor of 2 in both directions
if average_mesh && uniform
    [X, U] = compress_mesh(X, U);
end

% mean velocities and picture dimensions
mean_u = mean(u,3);
mean_v = mean(v,3);
dimensions = size(x);  

% Determine the number of images actually in memory
num_images = size(u,3);

% Determine resolution of velocity image
data_points = numel(x);

% determine units to be displayed in plots
if strcmp(xy_units, 'mm')
    x_dis = x*1000;
    y_dis = y*1000;
else
    x_dis = x;
    y_dis = y;
end

% Exactly define flow boundaries
[bnd_x, bnd_y, bnd_idx] = refine_bounds(x_dis, y_dis, u, v, mean_u, mean_v, direct, streamlines, open_flow, update_bnds);

% Filter raw images, to attempt to remove artifacts
if filter
    [u, v, mean_u, mean_v] = filter_images(u, v, bnd_idx, bnd_x, bnd_y);
end

% Find fluxating velocity of flow
flux_u = u - repmat(mean_u,1,1,num_images);
flux_v = v - repmat(mean_v,1,1,num_images);

% Remove excess portions of flow no longer in use
[mean_u, mean_v, flux_u, flux_v] = clip_bounds(bnd_idx, mean_u, mean_v, flux_u, flux_v);

% Calculate volume elements of the mesh
vol_frac = voln_piv_2D(x, y, bnd_idx, bnd_x, bnd_y);
    
clear u v

% Create a stacked data matrix for u and v velocities
flux_u      = reshape(flux_u, data_points, num_images);
flux_v      = reshape(flux_v, data_points, num_images);
mean_u      = reshape(mean_u, data_points, 1);
mean_v      = reshape(mean_v, data_points, 1);
vol_frac    = reshape(vol_frac, data_points, 1);

%% Perform Proper Orthogonal Decomposition

% Generate covariance matrix
covariance = cal_covariance_mat_2D(flux_u, flux_v, vol_frac, bnd_idx);

% Perform POD on fluctuating data
[pod_u, pod_v, lambda, modal_amp, cutoff] = POD_2D(covariance, flux_u, flux_v); 

% Figure out sign and apply flip
for i = 1:cutoff
    sign_flip = sign(mean(mean(sign(pod_v(:,i)./(pod_v(:,i) + eps)))));
    pod_u(:,i) = pod_u(:,i)*sign_flip;
    pod_v(:,i) = pod_v(:,i)*sign_flip;
    modal_amp(:,i) = modal_amp(:,i)*sign_flip;
end

close all

% Cluster resulting POD modes
if cluster
    [km_stoch, gm_stoch, gm_models, gm_groups, km_groups, centers] = ...
        cluster_POD(modal_amp, num_clusters, direct, save_figures);
end

% Calculate voritcity
[pod_vor, mean_vor] = calc_pod_vor(pod_u, pod_v, mean_u, mean_v, ...
                            dimensions, cutoff, bnd_idx, bnd_x, bnd_y, x, y);
                     
% Create a stacked data matrix for vorticity values
pod_vor  = reshape(pod_vor, data_points, cutoff);
mean_vor = reshape(mean_vor, data_points, 1);

%% Setup for Plotting and Plotting

% setup plot limit
if cutoff > 40
    num_plot = 40;
else
    num_plot = cutoff;
end

% Plot pod modes
if ~isempty(save_figures)
    
    % setup data structure
    data.x      = x_dis;
    data.y      = y_dis;
    data.bnd_idx = bnd_idx;
    
    % Plot vector field modes
    plot_vector_modes(data, pod_u, pod_v, num_plot, dimensions, 'POD', lambda, direct, save_figures, streamlines);
    
    % plot scalar field modes
    plot_scalar_modes(data, pod_u, num_plot, dimensions, 'u', lambda, direct, save_figures);
    plot_scalar_modes(data, pod_v, num_plot, dimensions, 'v', lambda, direct, save_figures);
    plot_scalar_modes(data, pod_vor, num_plot, dimensions, 'vorticity', lambda, direct, save_figures);
end

% Add mode zero
[modal_amp, lambda, pod_u, pod_v, pod_vor] = ...
    add_mode_zero(modal_amp, lambda, pod_u, pod_v, pod_vor, mean_u, mean_v, mean_vor);

%% Save / Return variables
run_num = floor(100000*rand(1));

% Place results in a structure
results_pod.run_num = run_num;
results_pod.exp_sampling_rate = exp_sampling_rate;

% mesh information
results_pod.x = x;
results_pod.y = y;
results_pod.x_dis = x_dis;
results_pod.y_dis = y_dis;
results_pod.dimensions = dimensions;
results_pod.bnd_idx = bnd_idx;
results_pod.bnd_x = bnd_x;
results_pod.bnd_y = bnd_y;
results_pod.uniform = uniform;
results_pod.vol_frac = vol_frac;

% fluctuating flow
results_pod.flux_u = flux_u;
results_pod.flux_v = flux_v;

% mean flow
results_pod.mean_u = mean_u;
results_pod.mean_v = mean_v;

% pod modes
results_pod.pod_u = pod_u;
results_pod.pod_v = pod_v;
results_pod.pod_vor = pod_vor;

% modal coordinates and energy captured, modes to 99%
results_pod.modal_amp = modal_amp;
results_pod.lambda = lambda;
results_pod.cutoff = cutoff;

% scaling factoring
results_pod.l_scale = l_scale;
results_pod.u_scale = u_scale;
results_pod.non_dim = non_dim;

if cluster
    % k-mean clustering data
    results_clust.km_groups = km_groups;
    results_clust.centers = centers;
    results_clust.km_stoch = km_stoch;
    
    % gaussian mixture model data
    results_clust.gm_models = gm_models;
    results_clust.gm_groups = gm_groups;
    results_clust.gm_stoch = gm_stoch;
    results_clust.cluster_range = 2:40;
end

% Save variables relavent to Galerkin to .mat files
if save_pod == true
    if cluster
        save([direct filesep 'POD Data' filesep 'POD_run_' num2str(run_num) '.mat'], ...
            'results_pod', 'results_clust', '-v7.3');
    else
        save([direct filesep 'POD Data' filesep 'POD_run_' num2str(run_num) '.mat'], ...
            'results_pod', '-v7.3');
    end
end

if nargout >= 1
    res_pod = results_pod;
end
if nargout == 2
    res_clust = results_clust;
end

% return format
format short g
end