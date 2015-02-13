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
% problem.num_modesG = 10
% Specify the number of modes that will calculated in Galerkin projection
%
% problem.plot_type = {'amp', 'fft'}
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
%
% problem.Re0_gen = @Re0_gen_shear
% Fuction handle used to produced Reynolds number for the run, takes one
% input of working project directory
%
% fft_window = [0 2000]
% Specify the hertz range that the fft plot should capture
%
% run_num = 'first'
% Specify which run this Galerkn Projection should be based from, default
% is to use the most recent


% Set format, clear figures, and set up correct directory
format long g
close all
clc;

% TODO need to allow to have multiple datasets present to more rapidly
% experiement with runs
% TODO make calculating coefficeints and plotting more generic so more can
% be added later

%List of fields that will be checked
fields = {  'num_modesG',     'plot_type',    'save_coef', ...
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
num_modesG      = problem.num_modesG;
run_num         = problem.run_num;
plot_type       = problem.plot_type;
save_coef       = problem.save_coef;
override_coef   = problem.override_coef;
tspan           = problem.tspan;
init            = problem.init;
direct          = problem.direct;
Re0_gen         = problem.Re0_gen;     
fft_window      = problem.fft_window;

% Check status of parrallel pool
if ~isempty(gcp('nocreate'))
    delete(gcp);
end

fprintf('\nLoading POD variables\n\n');

% Prompt User for folder if directory is not provided
if strcmp(direct, '');
    [direct_POD, direct] = prompt_folder('POD', run_num);
else
    [direct_POD, direct] = prompt_folder('POD', run_num, direct);
end

% Check folders are up to most recent format
update_folders(direct);

% Load POD variables
vars = load(direct_POD, 'results');

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
modal_amp_raw = vars.results.modal_amp_raw;   % modal amplitude of each image in pod basis
dimensions  = vars.results.dimensions;  % dimensions of mesh
vol_frac    = vars.results.vol_frac;    % mesh area size
bnd_idx     = vars.results.bnd_idx;     % location of boundaries
bnd_x       = vars.results.bnd_x;
bnd_y       = vars.results.bnd_y;
uniform     = vars.results.uniform;     % logical if mesh is uniform
run_num     = vars.results.run_num;     % POD run numbers
cutoff      = vars.results.cutoff;      % number of modes at cutoff

clear vars

% TODO may want to pass calculated high side velocity and calculate this
Re0 = Re0_gen(direct);      % Get Reynolds number        
z = ones(size(x));          % Depth of velocity field
t_scale = u_scale/l_scale;  % time scale
tspan = tspan*t_scale;      % non dimensionalized timescale

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
coef_problem.bnd_x          = bnd_x;
coef_problem.bnd_y          = bnd_y;
coef_problem.z              = z;
coef_problem.run_num        = run_num;
coef_problem.override_coef  = override_coef;
coef_problem.direct         = direct;
coef_problem.uniform        = uniform; %uniform;

% Generate unresolved coefficients
fprintf('Generating coefficients for unresolved modes using %d modes\n\n', cutoff);
[lc_dot , lc,  qc_2dot,  qc_dot,  qc]  = visocity_coefficients(coef_problem);
[lwc_dot, lwc, qwc_2dot, qwc_dot, qwc] = visocity_coefficients_ws(coef_problem); 

pod_ut = pod_u(:,1:num_modesG);
pod_vt = pod_v(:,1:num_modesG);

coef_problem.pod_u = pod_ut;
coef_problem.pod_v = pod_vt;

% Generate resolved coefficients
fprintf('Generating coefficients for resolved modes using %d modes\n\n', num_modesG);
[l_dot,  l,  q_2dot,  q_dot,  q] = visocity_coefficients(coef_problem);
[lw_dot, lw, qw_2dot, qw_dot, qw]= visocity_coefficients_ws(coef_problem);

% Free memory
clear pod_u pod_v mean_u mean_v

%% Modified coefficients 

fprintf('Calculating viscous disspation terms\n\n');

% Attempt to estimate the neglected viscoius dissapation function

%  calculate coefficeitns detailed by Noack
niu = viscious_dis(modal_amp_raw, num_modesG, lambda2, l, q_dot, q);
niuw = viscious_dis(modal_amp_raw, num_modesG, lambda2, lw, qw_dot, qw);

% calculate coefficients detailed by Couplet
niu_c = viscious_dis_couplet(modal_amp_raw, num_modesG, lc_dot, lc, qc_2dot, qc_dot, qc, Re0);
niuw_c = viscious_dis_couplet(modal_amp_raw, num_modesG, lwc_dot, lwc, qwc_2dot, qwc_dot, qwc, Re0);

                            
% From Couplet (v + v_tilde) ie (1/Re + v_tilde)
ni   = diag(niu) + (ones(num_modesG))/Re0;
ni_c = diag(niu_c) + (ones(num_modesG))/Re0; 
niw  = diag(niuw) + (ones(num_modesG))/Re0;
niw_c = diag(niuw_c) + (ones(num_modesG))/Re0;

% Calculate Constant terms terms
c    = l_dot/Re0 + q_2dot;
ci   = ni*l_dot + q_2dot;
ciw  = niw*lw_dot + qw_2dot;
ci_c = ni_c*l_dot+ q_2dot;
ciw_c = niw_c*lw_dot + qw_2dot;


% Calculate Linear terms
l    = l + q_dot;
li   = ni.*l   + q_dot;
liw  = niw*lw  + qw_dot;
li_c = ni_c.*l + q_dot;
liw_c = niw_c.*lw + qw_dot;

% System coefficients
Gal_coeff      = [c  l  q];
Gal_coeff_vis1 = [ci li q];
Gal_coeff_vis2 = [ci_c, li_c, q];
Gal_coeff_vis3 = [ciw, liw, q];
Gal_coeff_vis4 = [ciw_c, liw_c, q];

%% Time integration

% Will reduce number of coefficients if we want a smaller model
reduced_model_coeff_og   = -ode_coefficients(num_modesG, num_modesG, Gal_coeff);
reduced_model_coeff_vis1 = -ode_coefficients(num_modesG, num_modesG, Gal_coeff_vis1);
reduced_model_coeff_vis2 = -ode_coefficients(num_modesG, num_modesG, Gal_coeff_vis2);
reduced_model_coeff_vis3 = -ode_coefficients(num_modesG, num_modesG, Gal_coeff_vis3);
reduced_model_coeff_vis4 = -ode_coefficients(num_modesG, num_modesG, Gal_coeff_vis4);
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

fprintf('Performing ode113 on base Galerkin system\n');

% Integrate Base Galerkin System
tic;
[t1, modal_amp_og] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_og), tspan, ...
    modal_amp_raw(init,1:num_modesG), options);
toc1 = toc;
fprintf('Completed in %f6.4 seconds\n\n', toc1);



fprintf('Performing ode113 on Galerkin system with 1st viscous dissapation\n');

% Integrate Galerkin System with 1st viscous dissapation model
tic;
[t2, modal_amp_vis1] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_vis1), tspan, ...
    modal_amp_raw(init,1:num_modesG), options);
toc2 = toc;
fprintf('Completed in %f6.4 seconds\n\n', toc2);



fprintf('Performing ode113 on Galerkin system with 2nd viscous dissapation\n');

% Integrate Galerkin System with 2nd viscous dissapation model
tic;
[t3, modal_amp_vis2] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_vis2), tspan, ...
    modal_amp_raw(init,1:num_modesG), options);
toc3 = toc;
fprintf('Completed in %f6.4 seconds\n\n', toc3);




fprintf('Performing ode113 on Galerkin system with 2nd viscous dissapation\n');

% Integrate Galerkin System with 2nd viscous dissapation model
tic;
[t4, modal_amp_vis3] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_vis3), tspan, ...
    modal_amp_raw(init,1:num_modesG), options);
toc3 = toc;
fprintf('Completed in %f6.4 seconds\n\n', toc3);



fprintf('Performing ode113 on Galerkin system with 2nd viscous dissapation\n');

% Integrate Galerkin System with 2nd viscous dissapation model
tic;
[t5, modal_amp_vis4] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_vis4), tspan, ...
    modal_amp_raw(init,1:num_modesG), options);
toc3 = toc;
fprintf('Completed in %f6.4 seconds\n\n', toc3);


%% Plotting functions

% Prepare data
plot_data.num_modes     = num_modesG;
plot_data.direct        = direct;
plot_data.pod_ut        = pod_ut;
plot_data.pod_vt        = pod_vt;
plot_data.dimensions    = dimensions;
plot_data.fft_window    = fft_window;
plot_data.u_scale       = u_scale;
plot_data.l_scale       = l_scale;
plot_data.plot_type     = plot_type;
plot_data.sample_freq   = sample_freq;

all_ids = {'og', 'vis', 'vis2', 'vis3', 'vis4'};
all_modal_amps = {modal_amp_og, modal_amp_vis1, modal_amp_vis2, modal_amp_vis3, modal_amp_vis4};
all_t = {t1, t2, t3, t4, t5};

% Cycle through and plot all requested figures
for i = 1:size(all_ids,2);
   plot_data.id = all_ids{i};
   plot_data.modal_amp = all_modal_amps{i};
   plot_data.t = all_t{i};
   produce_plots(plot_data);
end

fprintf('Saving Galerkin Variables\n');

% Prepare data
results.num_run = run_num;
results.c = c;
results.ci = ci;
results.ci_c = ci_c;
results.l = l;
results.l_dot = l_dot;
results.li = li;
results.li_c = li_c;
results.ni = ni;
results.ni_c = ni_c;
results.q = q;
results.q_dot = q_dot;
results.q_2dot = q_2dot;
results.num_modesG = num_modesG;
results.modal_amp_og = modal_amp_og;
results.modal_amp_vis1 = modal_amp_vis1;
results.modal_amp_vis2 = modal_amp_vis2;
results.t1 = t1;
results.t2 = t2;
results.t3 = t3;
results.sample_freq = sample_freq;

% Save relavent coefficients
if save_coef == true
    save([direct '\Galerkin Coeff\Coeff_m' num2str(num_modesG) '_i' num2str(init) '_r'...
        num2str(run_num) '.mat'], 'results', '-v7.3');
end

if nargout == 1
    res = results;
end

% return format
format short g
end

