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
% problem.fft_window = [0 2000]
% Specify the hertz range that the fft plot should capture
%
% problem.run_num = 'first'
% Specify which run this Galerkn Projection should be based from, default
% is to use the most recent
%
% problem.dissapation = {'Noack', 'Couplet'}
% Specify which viscious dissapation method(s) to use
%
% problem.solution = {'Base', 'Weak'}
% Specify how visocity coefficients are calculated


% Set format, clear figures, and set up correct directory
format long g
close all
clc;

% TODO need to allow to have multiple datasets present to more rapidly
% experiement with runs
% TODO make calculating coefficeints and plotting more generic so more can
% be added later

%List of fields that will be checked
fields = {  'num_modesG',   'plot_type',    'save_coef', ...
            'override_coef','tspan',        'init', ...
            'direct' ,      'Re0_gen',      'fft_window', ...
            'run_num',      'dissapation',  'solution'};
        
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
solution        = problem.solution;
dissapation     = problem.dissapation;

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
bnd_x       = vars.results.bnd_x;       % location of flow boundaries normal to x
bnd_y       = vars.results.bnd_y;       % location of flow boundaries normal to y
uniform     = vars.results.uniform;     % logical if mesh is uniform
run_num     = vars.results.run_num;     % POD run numbers
cutoff      = vars.results.cutoff;      % number of modes at cutoff

clear vars

% TODO may want to pass calculated high side velocity and calculate this
Re0 = Re0_gen(direct, u_scale, l_scale);    % Get Reynolds number        
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
coef_problem.uniform        = uniform; 

% Prefill Cells
l       = cell(4,2);
l_dot   = cell(4,2);
q_2dot  = cell(4,2);
q_dot   = cell(4,2);
q       = cell(4,2);

% Generate unresolved coefficients, calculating only requested values
if ismember('Couplet', dissapation);
    if ismember('Base', solution)
        fprintf('Generating coefficients for unresolved modes using %d modes\n\n', cutoff);
        [l_dot{1,1} , l{1,1},  q_2dot{1,1},  q_dot{1,1},  q{1,1}]  = ...
            visocity_coefficients(coef_problem);
        l_dot{1,2}     = 'Base cutoff';
        l{1,2}         = 'Base cutoff';
        q_2dot{1,2}    = 'Base cutoff';
        q_dot{1,2}     = 'Base cutoff';
        q{1,2}         = 'Base cutoff';
    end
    if ismember('Weak', solution)
        fprintf('Generating coefficients for unresolved modes using %d modes for Weak Solution \n\n', cutoff);
        [l_dot{2,1} , l{2,1},  q_2dot{2,1},  q_dot{2,1},  q{2,1}] = ...
            visocity_coefficients_ws(coef_problem); 
        l_dot{2,2}     = 'Weak cutoff';
        l{2,2}         = 'Weak cutoff';
        q_2dot{2,2}    = 'Weak cutoff';
        q_dot{2,2}     = 'Weak cutoff';
        q{2,2}         = 'Weak cutoff';
    end
end

pod_ut = pod_u(:,1:num_modesG);
pod_vt = pod_v(:,1:num_modesG);

coef_problem.pod_u = pod_ut;
coef_problem.pod_v = pod_vt;

% Generate resolved coefficients, calculating only requested values
if ismember('Base', solution)
    fprintf('Generating coefficients for resolved modes using %d modes\n\n', num_modesG);
    [l_dot{3,1} , l{3,1},  q_2dot{3,1},  q_dot{3,1},  q{3,1}] = ... 
        visocity_coefficients(coef_problem);
    l_dot{3,2}     = 'Base num modes';
    l{3,2}         = 'Base num modes';
    q_2dot{3,2}    = 'Base num modes';
    q_dot{3,2}     = 'Base num modes';
    q{3,2}         = 'Base num modes';
end
if ismember('Weak', solution)
    fprintf('Generating coefficients for resolved modes using %d modes for Weak Solution \n\n', num_modesG);
    [l_dot{4,1} , l{4,1},  q_2dot{4,1},  q_dot{4,1},  q{4,1}] = ...
        visocity_coefficients_ws(coef_problem); 
    l_dot{4,2}     = 'Weak num modes';
    l{4,2}         = 'Weak num modes';
    q_2dot{4,2}    = 'Weak num modes';
    q_dot{4,2}     = 'Weak num modes';
    q{4,2}         = 'Weak num modes';
end

% Free memory
clear pod_u pod_v mean_u mean_v

%% Modified coefficients 

% Attempt to estimate the neglected viscoius dissapation function
fprintf('Calculating viscous disspation terms\n\n');

% Prefill cell
niu   = cell(4,2);

% calculate coefficients detailed by Couplet
if ismember('Couplet', dissapation);
    if ismember('Base', solution)
        niu{1,1} = viscious_dis_couplet(modal_amp_raw, num_modesG, l_dot{1,1}, ...
            l{1,1}, q_2dot{1,1}, q_dot{1,1}, q{1,1}, Re0);
        niu{1,2} = 'Base Couplet';
    end
    if ismember('Weak', solution)
        niu{2,1} = viscious_dis_couplet(modal_amp_raw, num_modesG, l_dot{2,1}, ...
            l{2,1}, q_2dot{2,1}, q_dot{2,1}, q{2,1}, Re0);
        niu{2,2} = 'Weak Couplet';
    end
end

%  calculate coefficeints detailed by Noack
if ismember('Base', solution)
    niu{3,1} = viscious_dis(modal_amp_raw, num_modesG, lambda2, l{3,1}, q_dot{3,1}, q{3,1});
    niu{3,2} = 'Base Noack';
end
if ismember('Weak', solution)
    niu{4,1} = viscious_dis(modal_amp_raw, num_modesG, lambda2, l{4,1}, q_dot{4,1}, q{4,1});
    niu{4,2} = 'Weak Noack';
end
                            
% From Couplet (v + v_tilde) ie (1/Re + v_tilde)
ni = cell(size(niu,1), 2);
for i = 1:size(niu,1);
    if ~isempty(niu{i,1})
        ni{i,1} = diag(niu{i,1}) + ones(num_modesG)/Re0;
        ni{i,2} = niu{i,2};
    end
end

% Calculate Constant terms terms
c_og   = l_dot{3,1}/Re0 + q_2dot{3,1};
ci  = cell(4,2);
for i = 1:size(ni,1);
    if ~isempty(ni{i,1})
        ci{i,1} = ni{i,1}*l_dot{i,1} + q_2dot{i,1};
        ci{i,2} = ni{i,2};
    end
end

% Calculate Linear terms
l_og    = l{3,1} + q_dot{3,1};
li = cell(4,2);
for i = 1:size(ni,1)
    if ~isempty(ni{i,1})
        li{i,1} = ni{i,1}*l{i,1} + q_dot{i,1};
        li{i,2} = ni{i,2};
    end
end

% System coefficients
q_og = q{3,1};
Gal_coeff      = [c_og  l_og  q_og];
Gal_coeff_vis = cell(4,2);
for i = 1:size(ni,1)
    if ~isempty(ni{i,1})
        Gal_coeff_vis{i,1} = [ci{i,1}, li{i,1}, q{i,1}];
        Gal_coeff_vis{i,2} = ni{i,2};
    end
end

%% Time integration

% Will reduce number of coefficients if we want a smaller model
reduced_model_coeff_og   = -ode_coefficients(num_modesG, num_modesG, Gal_coeff);
reduced_model_coeff_vis = cell(4,2);
for i = 1:size(ni,1)
    if ~isempty(Gal_coeff_vis{i,1})
        reduced_model_coeff_vis{i,1} = -ode_coefficients(num_modesG, num_modesG, Gal_coeff_vis{i,1});
        reduced_model_coeff_vis{i,2} = Gal_coeff_vis{i,2};
    end
end
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

% Integrate Base Galerkin System
fprintf('Performing ode113 on base Galerkin system\n');
tic;
[t_og, modal_amp_og] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_og), tspan, ...
    modal_amp_raw(init,1:num_modesG), options);
toc1 = toc;
fprintf('Completed in %f6.4 seconds\n\n', toc1);

% Integrate Galerkin System of requested methods
t_vis           = cell(4,2);
modal_amp_vis   = cell(4,2);
for i = 1:size(ni,1)
    if ~isempty(Gal_coeff_vis{i,1});
        fprintf('Performing ode113 on Galerkin system with %s\n', reduced_model_coeff_vis{i,2});
        tic;
        [t_vis{i,1}, modal_amp_vis{i,1}] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_vis{i,1}), tspan, ...
            modal_amp_raw(init,1:num_modesG), options);
        toc2 = toc;
        fprintf('Completed in %f6.4 seconds\n\n', toc2);
    end
end

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
plot_data.x             = x;
plot_data.y             = y;

all_ids = {'og'};
all_ids(2:5) = niu(:,2)';
all_modal_amps = {modal_amp_og};
all_modal_amps(2:5) =  modal_amp_vis(:,1)';
all_t = {t_og};
all_t(2:5) = t_vis(:,1)';

% Cycle through and plot all requested figures
for i = 1:size(all_ids,2);
    if ~isempty(all_ids{i})
        plot_data.id = all_ids{i};
        plot_data.modal_amp = all_modal_amps{i};
        plot_data.t = all_t{i};
        % TODO remove
        if i == 5
            plot_data.plot_type = {'amp', 'fft', 'video'};
        else
            plot_data.plot_type = {'amp', 'fft'};
        end
        produce_plots(plot_data);
    end
end

fprintf('Saving Galerkin Variables\n');

% Prepare data
results.num_run = run_num;
results.c_og = c_og;
results.l_og = l_og;
results.q_og = q_og;
results.l = l;
results.l_dot = l_dot;
results.q = q;
results.q_dot = q_dot;
results.q_2dot = q_2dot;
results.ci = ci;
results.li = li;
results.ni = ni;
results.niu = niu;
results.num_modesG = num_modesG;
results.modal_amp_og = modal_amp_og;
results.modal_amp_vis = modal_amp_vis;
results.t_og = t_og;
results.t_vis = t_vis;
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

