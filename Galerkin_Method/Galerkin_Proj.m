function res = Galerkin_Proj(varargin)
% GALERKIN_PROJ perform Galerkin projection of the POD system onto Navier
% Stokes. This requires POD_GEN to be run, Can produce save output graphs
% see below return results into RES if output is requested
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
% problem.use_chunks = false
% Set this to true if you are running out of memeory well write values of q
% to disk 



% Set format, clear figures, and set up correct directory
format long g
close all
clc;

% TODO need to allow to have multiple datasets present to more rapidly
% experiement with runs

%List of fields that will be checked
fields = {  'num_modesG',   'plot_type',    'save_coef', ...
            'override_coef','tspan',        'init', ...
            'direct' ,      'Re0_gen',      'fft_window', ...
            'run_num',      'dissapation',  'use_chunks'};
        
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
dissapation     = problem.dissapation;
use_chunks      = problem.use_chunks;

clear problem

% Check status of parrallel pool
if isempty(gcp('nocreate'));
    parpool('local');
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
pod_u       = vars.results.pod_u;       % streamwise pod modes
pod_v       = vars.results.pod_v;       % spanwise pod modes
pod_vor     = vars.results.pod_vor;     % vorticity modes
lambda      = vars.results.lambda;      % eigenvalues of modes
% cluster_range  = vars.results.cluster_range;    % number of variables in cluster
% centers     = vars.results.centers;     % cluster centers
dimensions  = vars.results.dimensions;  % dimensions of mesh
vol_frac    = vars.results.vol_frac;    % mesh area size
bnd_idx     = vars.results.bnd_idx;     % location of boundaries
bnd_x       = vars.results.bnd_x;       % location of flow boundaries normal to x
bnd_y       = vars.results.bnd_y;       % location of flow boundaries normal to y
uniform     = vars.results.uniform;     % logical if mesh is uniform
run_num     = vars.results.run_num;     % POD run numbers
cutoff      = vars.results.cutoff;    % number of modes at cutoff
modal_amp = vars.results.modal_amp;   % modal amplitude  from raw data


clear vars

% Get Reynolds number 
Re0 = Re0_gen(direct, u_scale, l_scale);      

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

% Ready Coef Problem Structure
coef_problem.x              = x;
coef_problem.y              = y;
coef_problem.use_chunks     = use_chunks;
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

% Free memory 
clear mean_u mean_v vol_frac bnd_x bnd_y 
clear z override_coef uniform

% Prefill Cells
lc = cell(2,2);
qc = cell(2,2);

% Generate values up to cutoff to be used for Couplet viscous dissipation
if ismember('Couplet', dissapation);
    % Calculate for traditional Galerkin Projection
    fprintf('Generating coefficients for unresolved modes using %d modes\n\n', cutoff);
    [lc{1,1}, qc{1,1}] = visocity_coefficients(coef_problem);
    lc{1,2}         = 'Base cutoff';
    qc{1,2}         = 'Base cutoff';
    
    % Calculate for weak formulation Galerkin Projection
    fprintf('Generating coefficients for unresolved modes using %d modes for Weak Solution \n\n', cutoff);
    lc{2,1}  = visocity_coefficients_ws(coef_problem); 
    lc{2,2}  = 'Weak cutoff';
    qc{2,1}  = qc{1,1}; 
    qc{2,2}  = 'Weak cutoff';
end

% determine number of models
base_models = 2;
linear_models = 8;
total_models = linear_models*2-base_models;

% Prefill cell
eddy= cell(total_models,2,length(num_modesG));
l   = cell(total_models,2,length(num_modesG));
q   = cell(total_models,2,length(num_modesG));
t           = cell(total_models,2,length(num_modesG));
modal_amp_sim   = cell(total_models,2,length(num_modesG));

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

for i = 1:length(num_modesG)
    close all
    
    % Pull current number of modes to be investigated
    num_modes = num_modesG(i);
    modal_TKE = sum(1/2*mean(modal_amp(:,2:num_modes).^2,1));
    
    % Created truncated pod basis
    pod_ut  = pod_u(:,1:num_modes);
    pod_vt  = pod_v(:,1:num_modes);
    pod_vort= pod_vor(:,1:num_modes);

    coef_problem.pod_u = pod_ut;
    coef_problem.pod_v = pod_vt;

    fprintf('Generating coefficients for resolved modes using %d modes\n\n', num_modes);
    [l{1,1,i}, q{1,1,i}] = visocity_coefficients(coef_problem);
    l{1,2,i} = 'Base Coeff';
    q{1,2,i} = 'Base Coeff';

    fprintf('Generating coefficients for resolved modes using %d modes for Weak Solution \n\n', num_modes);
    l{2,1,i} = visocity_coefficients_ws(coef_problem); 
    l{2,2,i} = 'Weak Coeff';
    q{2,1,i} = q{1,1,i};
    q{2,2,i} = 'Weak Coeff';
    
    % duplicate
    l(:,:,i) = repmat(l(1:2,:,i), total_models/2, 1, 1);
    q(:,:,i) = repmat(q(1:2,:,i), total_models/2, 1, 1);
    
%% Modified coefficients 

    % Attempt to estimate the neglected viscoius dissapation function
    fprintf('Calculating viscous dissapation terms\n\n');

    eddy(1,:,i) = {zeros(num_modes,1), 'Base Original'};
    eddy(2,:,i) = {zeros(num_modes,1), 'Weak Original'};
    
    % calculate coefficients detailed by Couplet
    if ismember('Couplet', dissapation);
        eddy{3,1,i} = viscious_dis_couplet(modal_amp, num_modes, lc{1,1}, qc{1,1}, Re0);
        eddy{3,2,i} = 'Modal Base Couplet';

        eddy{4,1,i} = viscious_dis_couplet(modal_amp, num_modes, lc{2,1}, qc{2,1}, Re0);
        eddy{4,2,i} = 'Modal Weak Couplet';
    end

    disp(eddy{3,1,i});
    
    % Calculate coefficeints detailed by Noack
    [eddy{5,1,i}, eddy{6,1,i}] = viscious_dis(modal_amp, num_modes, lambda, l{1,1,i}, q{1,1,i}, Re0);
    eddy{5,2,i} = 'Modal Base Noack';
    eddy{6,2,i} = 'Global Base Noack';

    [eddy{7,1,i}, eddy{8,1,i}] = viscious_dis(modal_amp, num_modes, lambda, l{2,1,i}, q{2,1,i}, Re0);
    eddy{7,2,i} = 'Modal Weak Noack';
    eddy{8,2,i} = 'Global Weak Noack';
    
    % duplicate TODO name this better
    eddy(linear_models+1:total_models,1,i) = eddy(base_models+1:linear_models,1,i);
    for j = base_models+1:linear_models
        eddy{j+linear_models-base_models,2,i} = ['NL ' eddy{j,2,i}];
    end
    
    % Perform final manipulation to prep integration
    [reduced_model_coeff] = integration_setup(eddy, Re0, l, q, i, total_models, linear_models, num_modes);

%% Time integration

    ao = modal_amp(init,1:num_modes);
    [t, modal_amp_sim] = time_integration(reduced_model_coeff, eddy, Re0, modal_TKE, ...
        i, t, modal_amp_sim, ao, tspan, total_models, linear_models, options);
    
%% Plotting functions

    % Prepare data
    plot_data.num_modes     = num_modes-1;
    plot_data.direct        = direct;
    plot_data.pod_ut        = pod_ut;
    plot_data.pod_vt        = pod_vt;
    plot_data.pod_vort      = pod_vort;
    plot_data.dimensions    = dimensions;
    plot_data.fft_window    = fft_window;
    plot_data.u_scale       = u_scale;
    plot_data.l_scale       = l_scale;
    plot_data.plot_type     = plot_type;
    plot_data.sample_freq   = sample_freq;
    plot_data.x             = x;
    plot_data.y             = y;
    plot_data.bnd_idx       = bnd_idx;

    all_t = t(:,1,i)';
    all_ids = eddy(:,2,i)';
    all_modal_amps =  modal_amp_sim(:,1,i)';

    % Cycle through and plot all requested figures
    for j = 1:total_models;
        if ~isempty(all_ids{j})
            plot_data.t = all_t{j};
            plot_data.id = all_ids{j};
            plot_data.modal_amp = all_modal_amps{j};
            produce_plots(plot_data);
        end
    end
end

fprintf('Saving Galerkin Variables\n');

% Prepare data
results.num_run = run_num;
results.l = l;
results.q = q;
results.eddy = eddy;
results.vis = 1/Re0;
results.num_modesG = num_modesG;
results.modal_amp_sim = modal_amp_sim;
results.t = t;
results.sample_freq = sample_freq;
results.linear_models = linear_models;
results.total_models = total_models;

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

