function Galerkin_Proj(varargin)
% GALERKIN_PROJ will calculate time coefficients for a provided POD basis
%
% GALERKIN_PROJ(NUM_PODS) generated time coefficients for NUM_PODS number
% of pod modes
%
% Galerkin_PROJ(NUM_PODS, PLOT_PRED) generated time coefficients for
% NUM_PODS number of pod modes. If plot_pred is true a short move is saved
% to the data folder of the predicted flow.

%%%%%%%%%%%%%%% VARIABLES USED FROM CHABOT MAIN CODE %%%%%%%%%%%%%%%%%%%%%%
%
% pod_u: pods in the u direction that have had thier sign flipped
% pod_v: pods in the v direction that have had their sign flipped
% mean_u: mean velocity in the u direction
% mean_v: mean velocity in the v direction
% dimensions: size of the piv images
% x: matrix of x values
% y: matrix of y values
% bnd_idx: matrix of -1's , 0's, 1's defining the detected boundaries

% Set format, clear figures, and set up correct directory
format long g
close all

%%%%%%%%%%%% ALTER THIS TO MAKE PORTABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd('D:\shear layer');

% set up function 
switch nargin
    case 0 
        % Default: run simulation for 10 pod modes, don't plot, run for
        % 100s, save coefficients
        num_pods  = 10;
        plot_pred = 'none';
        save_coef = true;
        tspan     = 0:0.01:100;
        init      = 1;
        direct    = '';
    case 1
        num_pods  = varargin{1};
    	plot_pred = 'none';
        save_coef = true;
        tspan     = 0:0.01:100;
        init      = 1;
        direct    = '';
    case 2
        num_pods  = varargin{1};
        plot_pred = varargin{2};
        save_coef = true;        
        tspan     = 0:0.01:100;
        init      = 1;
        direct    = '';  
    case 3 
        num_pods  = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        tspan     = 0:0.01:100;
        init      = 1;
        direct    = '';
    case 4
        num_pods  = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        tspan     = varargin{4};  
        init      = 1;
        direct    = '';
    case 5
        num_pods  = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        tspan     = varargin{4};  
        init      = varargin{5};
        direct    = '';
    case 6
        num_pods  = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        tspan     = varargin{4};  
        init      = varargin{5};
        direct    = varargin{6};
    otherwise
        error('Too many input arguments');
end

% Check that parallel pool is ready
if isempty(gcp)
    parpool;
end

gcp();

% Prompt User for folder if directory is not provided
if strcmp(direct, '');
    [data, direct] = prompt_folder('POD');
else
    [data, direct] = prompt_folder('POD', direct);
end
update_folders(direct);
load(data{1});


% TODO may want to pass calculated high side velocity and calculate this
Re0=0.28e6;         %Reynolds number based on separator plate length
z=ones(size(x));    %Depth of velocity field

% Determine sampling frequency from provided tspan
if length(tspan) > 2
    sample_freq = 1/(tspan(2) - tspan(1));
    disp(sample_freq);
else
    error('must provide tspan with a range');
end

% Subset of total pod modes for model generation
pod_ut = pod_u(:,1:num_pods);
pod_vt = pod_v(:,1:num_pods);

[lc_dot2, lc2, qc_2dot2, qc_dot2, qc2] = visocity_coefficients_fast(mean_u, mean_v, ...
    x, y, pod_u, pod_v, dimensions, vol_frac);

% coefficeints for the unresolved and resolved modes
[lc_dot, lc, qc_2dot, qc_dot, qc] = visocity_coefficients(mean_u, mean_v, ...
    x, y, pod_u, pod_v, dimensions, vol_frac, bnd_idx, z);

% Free memory
clear pod_u pod_v

% coefficients for the resolved modes
[l_dot, l, q_2dot, q_dot, q] = visocity_coefficients(mean_u, mean_v, ...
    x, y, pod_ut, pod_vt, dimensions, vol_frac, bnd_idx, z);

% Look more at this, 
niu = viscious_dis(eig_func, num_pods, lambda2, l, q_dot, q);
niu_c = viscious_dis_couplet(eig_func, num_pods, lambda2, ...
                                lc_dot, lc, qc_2dot, qc_dot, qc, Re0);

% From Couplet (v + v_tilde) ie (1/Re + v_tilde)
ni   = diag(niu) + (ones(num_pods))/Re0;
ni_c = diag(niu_c) + (ones(num_pods))/Re0; 

% Calculate Constant terms terms
c    = l_dot/Re0 + q_2dot;
ci   = ni*l_dot + q_2dot;
ci_c = ni_c*l_dot+ q_2dot;

% Calculate Linear terms
l    = l + q_dot;
li   = ni.*l   + q_dot;
li_c = ni_c.*l + q_dot;

% System coefficients
Gal_coeff      = [c  l  q];
Gal_coeff_vis1 = [ci li q];
Gal_coeff_vis2 = [ci_c, li_c, q];


% Will reduce number of coefficients if we want a smaller model
% reduced_model_coeff = -ode_coefficients(num_pods, num_pods, Gal_coeff);
reduced_model_coeff_vis1 = -ode_coefficients(num_pods, num_pods, Gal_coeff_vis1);
reduced_model_coeff_vis2 = -ode_coefficients(num_pods, num_pods, Gal_coeff_vis2);
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

% Integrate Galerkin System for no viscious damping, and the two viscous
% damping methods
% tic1 = tic;
% [t, modal_amp] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff), tspan, ...
%     eig_func(init,1:num_pods), options);
% toc(tic1);

tic2 = tic;
[t2, modal_amp_vis1] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_vis1), tspan, ...
    eig_func(init,1:num_pods), options);
toc(tic2);

tic3 = tic;
[t3, modal_amp_vis2] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_vis2), tspan, ...
    eig_func(init,1:num_pods), options);
toc(tic3);

%% Plotting functions

% Plot modal amplitudes
if any(strcmp(plot_pred, 'amp'))
%     plot_amp(modal_amp(:, 1:num_pods), t, direct, init);
    plot_amp(modal_amp_vis1(:, 1:num_pods), t2, direct, init, 'vis');
    plot_amp(modal_amp_vis2(:, 1:num_pods), t3, direct, init, 'vis');
end

% Produce time response video
if any(strcmp(plot_pred, 'video'))
%     plot_prediction(pod_ut, pod_vt, x, y, modal_amp, t, num_pods, dimensions, direct)
    plot_prediction(pod_ut, pod_vt, x, y, modal_amp_vis1, t2, num_pods, dimensions, direct)
    plot_prediction(pod_ut, pod_vt, x, y, modal_amp_vis2, t3, num_pods, dimensions, direct)
end

% Plot modal fft
if any(strcmp(plot_pred, 'fft'))
    if num_pods > 4
        num2plot = 1:4;
    else
        num2plot = 1:num_pods;
    end
    if size(t2,1) > 4096
        window_size = 8192;
    else
        window_size = size(t2,1);
    end
%     modal_fft(modal_amp, num2plot, size(pod_ut, 1), window_size, ...
%         sample_freq, [0 2000], direct);
    modal_fft(modal_amp_vis1, num2plot, size(pod_ut, 1), window_size,...
        sample_freq, [0 2000], direct, 'vis1')
    modal_fft(modal_amp_vis2, num2plot, size(pod_ut, 1), window_size,...
        sample_freq, [0 2000], direct, 'vis2')
end

% Save relavent coefficients
if save_coef == true
    save([direct '\Galerkin Coeff\Coeff_m' num2str(num_pods) 'i' num2str(init) '.mat'],...
        'ci', 'li',  'num_pods', 'modal_amp_vis1', 't2', 'l_dot', ...
        'q_2dot', 'q_dot', 'c', 'l', 'q', 'sample_freq');  % 't' ,'modal_amp', 
end

% return format
format short g


end

