function res = Galerkin_Proj(varargin)
% GALERKIN_PROJ perform Galerkin projection of the POD system onto Navier
% Stokes. This requires POD_GEN to be run, Can produce save output graphs
% see below
%
% GALERKIN_PROJ() prompt user for analysis folders for a given test run,
% will use all defaults detailed below
% 
% GALERKIN_PROJ(problem) Using fields provided in the structure PROBLEM
% sets up analysis specified by PROBLEM. all unfilled fields go to defaults
%
% problem.num_pods = 10
% Specify the number of modes that will calculated in Galerkin projection
%
% problem.plot_pred = {'amp', 'fft'}
% Specify which outout graphes are desired, current options are modal
% amplitude 'amp', fourier fast transform 'fft', and 'video which produces
% a video of the simulated flow
% 
% problem.save_coef = true
% Save relvant values to at .mat file
%
% problem.override_coef = false
% Overwrite previous coefficients for the same run number, if the run
% numbers are different new coefficients will be calculated
%
% problem.tspan = 0:0.01:100
% Specify time span for time integration
% 
% problem.init
% Speicify which image will constitute the initial conditions
%
% problem.direct = ''
% Specify directory that will be searched for POD data, default is to
% prompt user

% Set format, clear figures, and set up correct directory
format long g
close all
clc;

% TODO need to allow to have multiple datasets present to more rapidly
% experiement with runs
% TODO make calculating coefficeints and plotting more generic so more can
% be added later

%List of fields that will be checked
fields = {  'num_pods',     'plot_pred',    'save_coef', ...
            'override_coef','tspan',        'init', ...
            'direct' ,      'Re0_gen',      'fft_window', ...
            'run_num'};
        
% Parse problem structure provided to set it up correctly
if nargin == 1
    problem = parse_inputs(fields, @setdefaults_proj, varargin{1});
else
    fprintf('Provide a single structure as input, use help Galerkin_Proj for information.\n');
    fprintf('Using Defaults\n\n');
    problem = parse_inputs(fields, @setdefaults_proj);
end

% Create more readable names
num_pods        = problem.num_pods;
run_num         = problem.run_num;
plot_pred       = problem.plot_pred;
save_coef       = problem.save_coef;
override_coef   = problem.override_coef;
tspan           = problem.tspan;
init            = problem.init;
direct          = problem.direct;
Re0_gen         = problem.Re0_gen;     
fft_window      = problem.fft_window;

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
    [data, direct] = prompt_folder('POD', direct, run_num);
end

% Check folders are up to most recent format
update_folders(direct);

% Load POD variables
vars = load(data{1}, 'results');

% Create more readable names
x           = vars.results.x;           % mesh coordinates in x direction
y           = vars.results.y;           % mesh coordinates in y direction
u_scale     = vars.results.u_scale;     % velocity scaling
l_scale     = vars.results.l_scale;     % length scaling
mean_u      = vars.results.mean_u;      % mean streamwise velocity
mean_v      = vars.results.mean_v;      % mean spanwise velocity
pod_u       = vars.results.pod_u;       % streamwise pod modes
pod_v       = vars.results.pod_v;       % spanwise pod modes
lambda2     = vars.results.lambda2;     % eigenvalues of modes
modal_amp   = vars.results.modal_amp;    % Need to find a better name
dimensions  = vars.results.dimensions;  % dimensions of mesh
vol_frac    = vars.results.vol_frac;    % mesh area size
bnd_idx     = vars.results.bnd_idx;     % location of boundaries
uniform     = vars.results.uniform;     % logical if mesh is uniform
run_num     = vars.results.run_num;     % POD run numbers
cutoff      = vars.results.cutoff;      % number of modes at cutoff

clear results

% TODO may want to pass calculated high side velocity and calculate this
Re0 = Re0_gen();        % Get Reynolds number        
z = ones(size(x));      % Depth of velocity field

% Determine sampling frequency from provided tspan
if length(tspan) > 2
    sample_freq = 1/(tspan(2) - tspan(1));
    fprintf('Detected Sampling Frequency %6.2f Hz\n\n', sample_freq);
else
    error('must provide tspan with a range');
end

%% Calculate Coefficients
coef_problem.mean_u         = mean_u;
coef_problem.mean_v         = mean_v;
coef_problem.x              = x;
coef_problem.y              = y;
coef_problem.pod_u          = pod_u;
coef_problem.pod_v          = pod_v;
coef_problem.dimensions     = dimensions;
coef_problem.vol_frac       = vol_frac;
coef_problem.bnd_idx        = bnd_idx;
coef_problem.z              = z;
coef_problem.run_num        = run_num;
coef_problem.override_coef  = override_coef;
coef_problem.direct         = direct;
coef_problem.uniform        = uniform;

fprintf('Generating coefficients for unresolved modes using %d modes\n\n', cutoff);

[lc_dot, lc, qc_2dot, qc_dot, qc] = visocity_coefficients_new(coef_problem);

pod_ut = pod_u(:,1:num_pods);
pod_vt = pod_v(:,1:num_pods);

l_dot   = lc_dot(1:num_pods);
l       = lc(1:num_pods, 1:num_pods);
q_2dot  = qc_2dot(1:num_pods);
q_dot   = qc_dot(1:num_pods, 1:num_pods);
q       = qc(1:num_pods, 1:num_pods^2);

% Free memory
clear pod_u pod_v mean_u mean_v

%% Modified coefficients 

fprintf('Calculating viscous disspation terms\n\n');

% Attempt to estimate the neglected viscoius dissapation function

%  calculate coefficeitns detailed by Noack
niu = viscious_dis(modal_amp, num_pods, lambda2, l, q_dot, q);

% calculate coefficients detailed by Couplet
niu_c = viscious_dis_couplet(modal_amp, num_pods, lc_dot, lc, qc_2dot, qc_dot, qc, Re0, niu);

                            
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
    modal_amp(init,1:num_pods), options);
toc1 = toc;
fprintf('Completed in %f6.4 seconds\n\n', toc1);



fprintf('Performing ode113 on Galerkin system with 1st viscous dissapation\n');

% Integrate Galerkin System with 1st viscous dissapation model
tic;
[t2, modal_amp_vis1] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_vis1), tspan, ...
    modal_amp(init,1:num_pods), options);
toc2 = toc;
fprintf('Completed in %f6.4 seconds\n\n', toc2);



fprintf('Performing ode113 on Galerkin system with 2nd viscous dissapation\n');

% Integrate Galerkin System with 2nd viscous dissapation model
tic;
[t3, modal_amp_vis2] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_vis2), tspan, ...
    modal_amp(init,1:num_pods), options);
toc3 = toc;
fprintf('Completed in %f6.4 seconds\n\n', toc3);


%% Plotting functions

% Plot modal amplitudes
if any(strcmp(plot_pred, 'amp'))
    plot_amp(modal_amp(:, 1:num_pods), t1, direct, init);
    plot_amp(modal_amp_vis1(:, 1:num_pods), t2, direct, init, 'vis1');
    plot_amp(modal_amp_vis2(:, 1:num_pods), t3, direct, init, 'vis2');
end

% TODO significant overhaul to this function
% Produce time response video
if any(strcmp(plot_pred, 'video'))
    plot_prediction(pod_ut, pod_vt, x, y, modal_amp, t1, num_pods, dimensions, direct)
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
    modal_fft(modal_amp, num2plot, window_size(1), ...
        sample_freq, fft_window, direct);
    modal_fft(modal_amp_vis1, num2plot, window_size(2),...
        sample_freq, fft_window, direct, 'vis1')
    modal_fft(modal_amp_vis2, num2plot, window_size(3),...
        sample_freq, fft_window, direct, 'vis2')
end

fprintf('Saving Galerkin Variables\n');

results.num_run = run_num;
results.c = c;
results.ci = ci;
results.ci_c = ci_c;
results.l = l;
results.l_dot = l_dot;
results.li = li;
results.li_c = li_c;
results.q = q;
results.q_dot = q_dot;
results.q_2dot = q_2dot;
results.num_pods = num_pods;
results.modal_amp = modal_amp;
results.modal_amp_vis1 = modal_amp_vis1;
results.modal_amp_vis2 = modal_amp_vis2;
results.t1 = t1;
results.t2 = t2;
results.t3 = t3;
results.sample_freq = sample_freq;

% Save relavent coefficients
if save_coef == true
    save([direct '\Galerkin Coeff\Coeff_m' num2str(num_pods) '_i' num2str(init) '_r'...
        num2str(run_num) '.mat'], 'results', '-v7.3');
end

if nargout == 1
    res = results;
end

% return format
format short g
end

