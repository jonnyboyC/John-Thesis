function Galerkin_Proj(varargin)
% GALERKIN_PROJ will calculate time coefficients for a provided POD basis
%
% GALERKIN_PROJ(NUM_PODS) generated time coefficients for NUM_PODS number
% of pod modes
%
% Galerkin_PROJ(NUM_PODS, PLOT_PRED) generated time coefficients for
% NUM_PODS number of pod modes. If plot_pred is true a short move is saved
% to the data folder of the predicted flow.


% Set format, clear figures, and set up correct directory
format long g
close all
clc;

%TODO major overhaul, look to make structure set up function 
switch nargin
    case 0 
        % Default: run simulation for 10 pod modes, don't plot, run for
        % 100s, save coefficients
        num_pods  = 10;
        plot_pred = 'none';
        save_coef = true;
        over_coef = false;
        tspan     = 0:0.01:100;
        init      = 1;
        direct    = '';
    case 1
        num_pods  = varargin{1};
    	plot_pred = 'none';
        save_coef = true;
        over_coef = false;
        tspan     = 0:0.01:100;
        init      = 1;
        direct    = '';
    case 2
        num_pods  = varargin{1};
        plot_pred = varargin{2};
        save_coef = true;      
        over_coef = false;
        tspan     = 0:0.01:100;
        init      = 1;
        direct    = '';  
    case 3 
        num_pods  = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        over_coef = false;
        tspan     = 0:0.01:100;
        init      = 1;
        direct    = '';
    case 4
        num_pods  = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        over_coef = varargin{4};
        tspan     = 0:0.01:100;  
        init      = 1;
        direct    = '';
    case 5
        num_pods  = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        over_coef = varargin{4};
        tspan     = varargin{5}; 
        init      = 1;
        direct    = '';
    case 6
        num_pods  = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        over_coef = varargin{4};
        tspan     = varargin{5};  
        init      = varargin{6};
        direct    = '';
    case 7
        num_pods  = varargin{1};
        plot_pred = varargin{2};
        save_coef = varargin{3};
        over_coef = varargin{4};
        tspan     = varargin{5};  
        init      = varargin{6};
        direct    = varargin{7};
    otherwise
        error('Too many input arguments');
end

% Check status of parrallel pool
if isempty(gcp)
    parpool;
end

gcp();

fprintf('\nLoading POD variables\n\n');

% Prompt User for folder if directory is not provided
if strcmp(direct, '');
    [data, direct] = prompt_folder('POD');
else
    [data, direct] = prompt_folder('POD', direct);
end

% Check folders are up to most recent format
update_folders(direct);

% Load POD variables
vars = load(data{1}, 'x', 'y', 'u_scale', 'l_scale', 'mean_u', 'mean_v', 'pod_u', ...
    'pod_v', 'eig_func', 'lambda2', 'dimensions', 'vol_frac', 'bnd_idx', 'uniform', ...
    'run_num', 'cutoff');

x           = vars.x;           % mesh coordinates in x direction
y           = vars.y;           % mesh coordinates in y direction
u_scale     = vars.u_scale;     % velocity scaling
l_scale     = vars.l_scale;     % length scaling
mean_u      = vars.mean_u;      % mean streamwise velocity
mean_v      = vars.mean_v;      % mean spanwise velocity
pod_u       = vars.pod_u;       % streamwise pod modes
pod_v       = vars.pod_v;       % spanwise pod modes
lambda2     = vars.lambda2;     % eigenvalues of modes
eig_func    = vars.eig_func;    % Need to find a better name
dimensions  = vars.dimensions;  % dimensions of mesh
vol_frac    = vars.vol_frac;    % mesh area size
bnd_idx     = vars.bnd_idx;     % location of boundaries
uniform     = vars.uniform;     % logical if mesh is uniform
run_num     = vars.run_num;     % run number of main_code_chabot
cutoff      = vars.cutoff;      % number of modes at cutoff


% TODO may want to pass calculated high side velocity and calculate this
Re0=0.28e6;                 % Reynolds number based on separator plate length   
z=ones(size(x));            % Depth of velocity field
% v = abs((u_scale*l_scale))/Re0;  % Vicosity

% Determine sampling frequency from provided tspan
if length(tspan) > 2
    sample_freq = 1/(tspan(2) - tspan(1));
    fprintf('Detected Sampling Frequency %6.2f Hz\n\n', sample_freq);
else
    error('must provide tspan with a range');
end

% Subset of total pod modes for model generation
pod_ut = pod_u(:,1:num_pods);
pod_vt = pod_v(:,1:num_pods);

%% Base Coefficients

fprintf('Generating coefficients for unresolved modes using %d modes\n\n', cutoff);

% coefficeints for the unresolved and resolved modes
if uniform == true
    % Use build in laplcian, and gradient functions
    [lc_dot, lc, qc_2dot, qc_dot, qc] = visocity_coefficients_fast(mean_u, mean_v, ...
        x, y, pod_u, pod_v, dimensions, vol_frac, run_num, over_coef, direct);
else
    % Old method allows for non_uniform mesh, SLOW
    [lc_dot, lc, qc_2dot, qc_dot, qc] = visocity_coefficients(mean_u, mean_v, ...
        x, y, pod_u, pod_v, dimensions, vol_frac, bnd_idx, z, run_num, over_coef, direct);
end

% Free memory
clear pod_u pod_v

fprintf('Generating coefficients for resolved modes using %d modes\n\n', num_pods);

% coefficients for the resolved modes
if uniform == true
    % Use build in laplcian, and gradient functions
    [l_dot, l, q_2dot, q_dot, q] = visocity_coefficients_fast(mean_u, mean_v, ...
        x, y, pod_ut, pod_vt, dimensions, vol_frac);
else
    % Old method allows for non_uniform mesh, SLOW
    [l_dot, l, q_2dot, q_dot, q] = visocity_coefficients(mean_u, mean_v, ...
        x, y, pod_ut, pod_vt, dimensions, vol_frac, bnd_idx, z);
end

%% Modified coefficients 

fprintf('Calculating viscous disspation terms\n\n');

% Attempt to estimate the neglected viscoius dissapation function
niu = viscious_dis(eig_func, num_pods, lambda2, l, q_dot, q);
niu_c = viscious_dis_couplet(eig_func, num_pods, ...
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

%% Time integration

% Will reduce number of coefficients if we want a smaller model
reduced_model_coeff = -ode_coefficients(num_pods, num_pods, Gal_coeff);
reduced_model_coeff_vis1 = -ode_coefficients(num_pods, num_pods, Gal_coeff_vis1);
reduced_model_coeff_vis2 = -ode_coefficients(num_pods, num_pods, Gal_coeff_vis2);
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

fprintf('Performing ode113 on base Galerkin system\n');

% Integrate Base Galerkin System
tic;
[t1, modal_amp] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff), tspan, ...
    eig_func(init,1:num_pods), options);
toc1 = toc;
fprintf('Completed in %f6.4 seconds\n\n', toc1);



fprintf('Performing ode113 on Galerkin system with 1st viscous dissapation\n');

% Integrate Galerkin System with 1st viscous dissapation model
tic;
[t2, modal_amp_vis1] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_vis1), tspan, ...
    eig_func(init,1:num_pods), options);
toc2 = toc;
fprintf('Completed in %f6.4 seconds\n\n', toc2);



fprintf('Performing ode113 on Galerkin system with 2nd viscous dissapation\n');

% Integrate Galerkin System with 2nd viscous dissapation model
tic;
[t3, modal_amp_vis2] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_vis2), tspan, ...
    eig_func(init,1:num_pods), options);
toc3 = toc;
fprintf('Completed in %f6.4 seconds\n\n', toc3);


%% Plotting functions

% Plot modal amplitudes
if any(strcmp(plot_pred, 'amp'))
    plot_amp(modal_amp(:, 1:num_pods), t1, direct, init);
    plot_amp(modal_amp_vis1(:, 1:num_pods), t2, direct, init, 'vis1');
    plot_amp(modal_amp_vis2(:, 1:num_pods), t3, direct, init, 'vis2');
end

% Produce time response video
if any(strcmp(plot_pred, 'video'))
    plot_prediction(pod_ut, pod_vt, x, y, modal_amp, t, num_pods, dimensions, direct)
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
    window_size = zeros(3,1);
    t = {t1, t2, t3};
    sample_freq = sample_freq*(u_scale/l_scale);
    for i = 1:3
        t{i} = t{i}*(u_scale/l_scale);
        if size(t{i},1) > 4096
            window_size(i) = 20000; %8192
        else
            window_size(i) = size(t{i},1);
        end
    end
%     modal_fft(modal_amp, num2plot, window_size(1), ...
%         sample_freq, [0 2000], direct);
    modal_fft(modal_amp_vis1, num2plot, window_size(2),...
        sample_freq, [0 2000], direct, 'vis1')
%     modal_fft(modal_amp_vis2, num2plot, window_size(3),...
%         sample_freq, [0 2000], direct, 'vis2')
end

fprintf('Saving Galerkin Variables\n');

% Save relavent coefficients
if save_coef == true
    save([direct '\Galerkin Coeff\Coeff_m' num2str(num_pods) 'i' num2str(init) '.mat'],...
        'c', 'ci', 'ci_c', 'l', 'l_dot', 'li', 'li_c', 'q', 'num_pods', 'q_2dot', 'q_dot',...
        'modal_amp', 'modal_amp_vis1', 'modal_amp_vis2', 't1', 't2', 't3', 'sample_freq', ...
        '-v7.3');  % 
end

% return format
format short g


end

