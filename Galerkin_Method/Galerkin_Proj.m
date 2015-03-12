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
pod_vor     = vars.results.pod_vor;     % vorticity modes
%lambda2     = vars.results.lambda2;     % eigenvalues of modes
modal_amp_mean = vars.results.modal_amp_mean;   % modal amplitude mean from raw data
modal_amp_flux = vars.results.modal_amp_flux;   % modal amplitude flucuations from raw data
dimensions  = vars.results.dimensions;  % dimensions of mesh
vol_frac    = vars.results.vol_frac;    % mesh area size
bnd_idx     = vars.results.bnd_idx;     % location of boundaries
bnd_x       = vars.results.bnd_x;       % location of flow boundaries normal to x
bnd_y       = vars.results.bnd_y;       % location of flow boundaries normal to y
uniform     = vars.results.uniform;     % logical if mesh is uniform
run_num     = vars.results.run_num;     % POD run numbers
cutoff      = vars.results.cutoff;    % number of modes at cutoff

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

% Add mode 0
pod_u = [mean_u, pod_u];
pod_v = [mean_v, pod_v];

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
num_models = 6;

% Prefill cell
niu = cell(num_models,2,length(num_modesG));
ni  = cell(num_models,2,length(num_modesG));
l   = cell(num_models,2,length(num_modesG));
li  = cell(num_models,2,length(num_modesG));
q   = cell(num_models,2,length(num_modesG));
t           = cell(num_models,2,length(num_modesG));
Gal_coeff   = cell(num_models,2,length(num_modesG));
modal_amp   = cell(num_models,2,length(num_modesG));

options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

for i = 1:length(num_modesG)
    close all
    
    % Pull current number of modes to be investigated
    num_modes = num_modesG(i);
    modal_TKE = sum(1/2*mean(modal_amp_flux(:,2:num_modes).^2,1));
    
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
    
    l(:,:,i) = repmat(l(1:2,:,i), num_models/2, 1, 1);
    q(:,:,i) = repmat(q(1:2,:,i), num_models/2, 1, 1);
    
%% Modified coefficients 

    % Attempt to estimate the neglected viscoius dissapation function
    fprintf('Calculating viscous disspation terms\n\n');

    niu(1,:,i) = {zeros(num_modes,1), 'Base Original'};
    niu(2,:,i) = {zeros(num_modes,1), 'Weak Original'};
    
    % calculate coefficients detailed by Couplet
    if ismember('Couplet', dissapation);
        niu{3,1,i} = viscious_dis_couplet(modal_amp_flux, modal_amp_mean, num_modes, lc{1,1}, qc{1,1}, Re0);
        niu{3,2,i} = 'Base Couplet';

        niu{4,1,i} = viscious_dis_couplet(modal_amp_flux, modal_amp_mean, num_modes, lc{2,1}, qc{2,1}, Re0);
        niu{4,2,i} = 'Weak Couplet';
    end

    % Calculate coefficeints detailed by Noack
    niu{5,1,i} = viscious_dis(modal_amp_flux, modal_amp_mean, num_modes, l{1,1,i}, q{1,1,i}, Re0);
    niu{5,2,i} = 'Base Noack';

    niu{6,1,i} = viscious_dis(modal_amp_flux, modal_amp_mean, num_modes, l{2,1,i}, q{2,1,i}, Re0);
    niu{6,2,i} = 'Weak Noack';
    
    reduced_model_coeff_vis = cell(num_models,2);

    % Final setup for time integration
    for j = 1:num_models
        if ~isempty(niu{j,-6,i})
            % Calculate Total Visocity
            ni{j,1,i} = niu{j,1,i} + 1/Re0;
            ni{j,2,i} = niu{j,2,i};
            
            % Setup up system coefficients
            Gal_coeff{j,1,i} = [l{j,1,i}, q{j,1,i}];
            Gal_coeff{j,2,i} = l{j,2,i};

            % rearrage into 1D form
            reduced_model_coeff_vis{j,1} = ode_coefficients(num_modes, Gal_coeff{j,1,i});
            reduced_model_coeff_vis{j,2} = Gal_coeff{j,2,i};
        end
    end

%% Time integration
    % Integrate Galerkin System of requested methods
    for j = 1:size(ni,1)
        if ~isempty(Gal_coeff{j,1});
            fprintf('Performing ode113 on Galerkin system with %s\n', reduced_model_coeff_vis{j,2});
            ao = modal_amp_flux(init,1:num_modes)+modal_amp_mean(init, 1:num_modes);
            tic;
            [t{j,1,i}, modal_amp{j,1,i}] = ode113(@(t,y) system_odes2(t, y, reduced_model_coeff_vis{j,1}, ni{j,1,i}, modal_TKE, false), ...
                tspan, ao, options);
            toc2 = toc;

            fprintf('Completed in %f6.4 seconds\n\n', toc2);
            
            modal_amp{j,1,i} = modal_amp{j,1,i}(:,2:end);
            
            t{j,2,i} = reduced_model_coeff_vis{j,2};
            modal_amp{j,2,i}  = reduced_model_coeff_vis{j,2};
        end
    end
%% Plotting functions

    % Prepare data
    plot_data.num_modes     = num_modes-1;
    plot_data.direct        = direct;
    plot_data.pod_ut        = pod_ut;
    plot_data.pod_vt        = pod_vt;
    plot_data.pod_vor       = pod_vort;
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
    all_ids = niu(:,2,i)';
    all_modal_amps =  modal_amp(:,1,i)';

    % Cycle through and plot all requested figures
    for j = 1:num_models;
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
results.ni = ni;
results.niu = niu;
results.num_modesG = num_modesG;
results.modal_amp = modal_amp;
results.t = t;
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

