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
% pod_u1: pods in the u direction that have had thier sign flipped
% pod_v1: pods in the v direction that have had their sign flipped
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
        num_pods = 10;
        plot_pred = 'none';
        save_coef = true;
        tspan = [0:0.01:100];
        init = 1;
        direct = '';
    case 1
        num_pods = varargin{1};
    	plot_pred = 'none';
        save_coef = true;
        tspan = [0:0.01:100];
        init = 1;
        direct = '';
    case 2
        num_pods = varargin{1};
        plot_pred = varargin{2};
        save_coef = true;        
        tspan = [0:0.01:100];
        init = 1;
        direct = '';  
    case 3 
        num_pods = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        tspan = [0:0.01:100];
        init = 1;
        direct = '';
    case 4
        num_pods = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        tspan = varargin{4};  
        init = 1;
        direct = '';
    case 5
        num_pods = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        tspan = varargin{4};  
        init = varargin{5};
        direct = '';
    case 6
        num_pods = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        tspan = varargin{4};  
        init = varargin{5};
        direct = varargin{6};
    otherwise
        error('Too many input arguments');
end

% Check that parallel pool is ready
if isempty(gcp)
    parpool;
end

pool = gcp;

% Prompt User for folder if directory is not provided
if strcmp(direct, '');
    [data, direct] = prompt_folder('POD');
else
    [data, direct] = prompt_folder('POD', direct);
end
update_folders(direct);
load(data{1});

Re0=0.28e6;         %Reynolds number based on separator plate length
z=ones(size(x));    %Depth of velocity field

if length(tspan) > 2
    sample_freq = 1/(tspan(2) - tspan(1));
    disp(sample_freq);
else
    error('must provide tspan with a range');
end

% TODO if we decide we want to calculate for range of number of pod modes,
% place for loop here

% Take a truncation of calculated pods
pod_ut = pod_u1(:,1:num_pods);
pod_vt = pod_v1(:,1:num_pods);

[l_dot, l, q_2dot, q_dot, q] = visocity_coefficients(mean_u, mean_v, ...
    x, y, pod_ut, pod_vt, dimensions, vol_frac, bnd_idx, z);

% Look more at this, 
niu = viscious_dis(eig_func, num_pods, lambda2, l, q_dot, q);
ni  = diag(niu) + (ones(num_pods)-eye(num_pods))/Re0;

% c  =  l_dot/Re0 + q_2dot;
ci = l_dot/Re0.*niu+q_2dot;

li = ni.*l+q_dot;
% l  = l+q_dot;

% Gal_coeff     = [c  l  q];
Gal_coeff_vis = [ci li q];


% Will reduce number of coefficients if we want a smaller model
% reduced_model_coeff = -ode_coefficients(num_pods, num_pods, Gal_coeff);
reduced_model_coeff_vis = -ode_coefficients(num_pods, num_pods, Gal_coeff_vis);
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

% TODO investigate situations where various ode solvers are faster
% tic1 = tic;
% [t, modal_amp] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff), tspan, ...
%     eig_func(init,1:num_pods), options);
% toc(tic1);

tic2 = tic;
[t2, modal_amp_vis] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_vis), tspan, ...
    eig_func(init,1:num_pods), options);
toc(tic2);

% Provide only modal flucations ie only turblent portion of modes
% TODO check the validity of this statement
% modal_amp = modal_amp - ones(size(modal_amp,1), 1)*mean(modal_amp);

% TODO update for plot_prediction
if any(strcmp(plot_pred, 'amp'))
%     plot_amp(modal_amp(:, 1:num_pods), t, direct, init);
    plot_amp(modal_amp_vis(:, 1:num_pods), t2, direct, init, 'vis');
end
if any(strcmp(plot_pred, 'video'))
%     plot_prediction(pod_ut, pod_vt, x, y, modal_amp, t, num_pods, dimensions, direct)
    plot_prediction(pod_ut, pod_vt, x, y, modal_amp_vis, t2, num_pods, dimensions, direct)
end
if any(strcmp(plot_pred, 'fft'))
    if num_pods > 4
        num2plot = 1:4;
    else
        num2plot = 1:num_pods;
    end
    if size(t2,1) > 4096
        window_size = 4096;
    else
        window_size = size(t2,1);
    end
%     modal_fft(modal_amp, num2plot, size(pod_ut, 1), window_size, ...
%         sample_freq, [0 2000], direct);
    modal_fft(modal_amp_vis, num2plot, size(pod_ut, 1), window_size,...
        sample_freq, [0 2000], direct, 'vis')
end

if save_coef == true
    save([direct '\Galerkin Coeff\Coeff_m' num2str(num_pods) 'i' num2str(init) '.mat'],...
        'ci', 'li',  'num_pods', 'modal_amp_vis', 't2', 'l_dot', ...
        'q_2dot', 'q_dot', 'q', 'sample_freq');  % 't' ,'modal_amp', 'c', 'l',
end

