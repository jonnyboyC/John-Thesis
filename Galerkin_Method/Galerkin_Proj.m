function [res_coef, res_int] = Galerkin_Proj(varargin)
% GALERKIN_PROJ perform Galerkin projection of the POD system onto Navier
% Stokes. This requires POD_GEN to be run, can produce and save graphs
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
% problem.plot_type = {'amp', 'fft', 'energy'}
% Specify which outout graphes are desired, current options are modal
% amplitude 'amp', fourier fast transform 'fft', system energy 'energy', 
% 'video' which produces a video of the simulated flow using quivers, 
% 'video stream' will produce a video of simulated flow using streamlines,
% caution very slow
% 
% problem.save_coef = true
% Save relvant values to at .mat file
%
% problem.override_coef = false
% Overwrite previous coefficients for the same run number, if the run
% numbers are different new coefficients will be calcfinteulated
%
% problem.tspan = 0:0.0001:1
% Specify time span for time integration
% 
% problem.init = 1
% Specify which image will constitute the initial conditions
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
% problem.dissapation = {'Least Squares', 'Averaged'}
% Specify which viscious dissapation method(s) to use
%
% problem.time_int = true
% If false no time integration or plotting will be performed, used to
% generated the Galerkin systems for MOD_POD
%
% problem.calc_coef = true
% If false attempt to locate previous runs files to time integration
%
% problem.use_chunks = false
% Set this to true if you are running out of memeory well write values of q
% to disk 
%
% problem.classify_sim = true
% Classify simulated results to emprical data

% Set format, clear figures, and set up correct directory
format long g
close all
clc;

%List of fields that will be checked
fields = {  'num_modesG',   'plot_type',    'save_coef', ...
            'override_coef','tspan',        'init', ...
            'direct' ,      'Re0_gen',      'fft_window', ...
            'run_num',      'dissapation',  'time_int', ...
            'use_chunks',   'calc_coef',    'classify_sim'};
        
% Parse problem structure provided to set it up correctly
if nargin == 1
    problem = parse_inputs(fields, @setdefaults_proj, varargin{1});
else
    fprintf('Provide a single structure as input, use help Galerkin_Proj for information.\n');
    fprintf('Using Defaults\n\n');
    problem = parse_inputs(fields, @setdefaults_proj);
end

% Create more readable names
num_modesG      = problem.num_modesG; % Add mean flow
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
calc_coef       = problem.calc_coef;
time_int        = problem.time_int;
classify_sim    = problem.classify_sim;
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
vars = load(direct_POD, 'results_pod');

% Create more readable names
X           = vars.results_pod.X;           % mesh coordinates
u_scale     = vars.results_pod.u_scale;     % velocity scaling
l_scale     = vars.results_pod.l_scale;     % length scaling
pod_U       = vars.results_pod.pod_U;       % pod modes
pod_W       = vars.results_pod.pod_W;       % vorticity pod modes
lambda      = vars.results_pod.lambda;      % eigenvalues of modes
dimensions  = vars.results_pod.dimensions;  % dimensions of mesh
bnd_idx     = vars.results_pod.bnd_idx;     % location of boundaries
non_dim     = vars.results_pod.non_dim;     % was the result non_dimensionalized
run_num     = vars.results_pod.run_num;     % POD run numbers
cutoff      = vars.results_pod.cutoff;      % number of modes at cutoff
modal_amp   = vars.results_pod.modal_amp;   % modal amplitude  from raw data

if calc_coef
    uniform     = vars.results_pod.uniform;      % Is the mesh uniform 
    bnd_X       = vars.results_pod.bnd_X;       % location of flow boundaries normal to x
    vol_frac    = vars.results_pod.vol_frac;    % mesh area size
end

% Free memory
clear vars

if time_int && classify_sim
    % Load Cluster variables
    vars = load(direct_POD, 'results_clust');
    
    centers         = vars.results_clust.centers;     % k-means cluster centers
    km_stoch        = vars.results_clust.km_stoch;    % k-mean stochastic matrix
    gm_stoch        = vars.results_clust.gm_stoch;    % gaussian mixture stochastic matrix
    gm_models       = vars.results_clust.gm_models;   % gaussian mixture models
    cluster_range   = vars.results_clust.cluster_range;    % number of variables in cluster
    
    % Free memory
    clear vars
end


% Determine Reynolds number 
Re0 = Re0_gen(direct);  

% Determine kinematic visocity
vis = calc_vis(u_scale, l_scale, Re0, non_dim);

% Determine sampling frequency from provided tspan
if length(tspan) > 2
    sample_freq = 1/(tspan(2) - tspan(1));
    fprintf('Detected Sampling Frequency %6.2f Hz\n\n', sample_freq);
else
    error('must provide tspan with a range');
end

% Create modify time scale if non-dimenionalized
[t_scale, tspan] = calc_t_scale(u_scale, l_scale, non_dim, tspan);

% Preallocate for future save function
futures = 0;

%% Setup variables to Coefficient Calculations

if calc_coef 
    
    % Ready Coef Problem Structure
    coef_problem.X              = X;
    coef_problem.use_chunks     = use_chunks;
    coef_problem.pod_U          = pod_U;
    coef_problem.dimensions     = dimensions;
    coef_problem.vol_frac       = vol_frac;
    coef_problem.bnd_idx        = bnd_idx;
    coef_problem.bnd_X          = bnd_X;
    coef_problem.run_num        = run_num;
    coef_problem.override_coef  = override_coef;
    coef_problem.direct         = direct;
    coef_problem.custom         = false;
    coef_problem.uniform        = uniform;
    
    % Free memory
    clear vol_frac bnd_X 
    
    % Prefill Cells
    lc = cell(2,2);
    qc = cell(2,2);
    
    % Generate values up to cutoff to be used for Couplet viscous dissipation
    if ismember('Least Squares', dissapation);
        fprintf('Generating coefficients for unresolved modes using %d modes\n', cutoff);
        
        % Calculate for traditional Galerkin Projection
        [lc{1,1}, qc{1,1}] = galerkin_coefficients(coef_problem);
        lc{1,2}  = 'Base cutoff';
        qc{1,2}  = 'Base cutoff';
        
        % Calculate for weak formulation Galerkin Projection
        lc{2,1}  = galerkin_coefficients_ws(coef_problem);
        lc{2,2}  = 'Weak cutoff';
        qc{2,1}  = qc{1,1};
        qc{2,2}  = 'Weak cutoff';
    end
end

% TODO eventually make this more dynamic
% determine number of models
base_models     = 2;
linear_models   = 8;
total_models    = linear_models*2-base_models;

% Initialize coefficeints
eddy    = cell(total_models,2,length(num_modesG));
l       = cell(total_models,2,length(num_modesG));
q       = cell(total_models,2,length(num_modesG));

% Initialize if time integration is being performed
if time_int
    t               = cell(total_models,2,length(num_modesG));
    modal_amp_sim   = cell(total_models,2,length(num_modesG));
end

% Initialize scores if clustering
if classify_sim
    scores_km       = cell(total_models,2,length(num_modesG));
    scores_gm       = cell(total_models,2,length(num_modesG));
end

% Perform calculation on each set of modes requested
for i = 1:length(num_modesG)
    fprintf('\n');
    
    % Determine if custom modes have been used, and setup appropriately 
    if iscell(num_modesG)
        custom = true;
        num_modes = length(num_modesG{i})+1;
        modes = num_modesG{i}+1;
    else
        custom = false;
        num_modes = num_modesG(i)+1;
        modes = 2:num_modes;
    end
    
    % Calculate empirical average TKE for selected modes
    modal_TKE = mean(sum(1/2*modal_amp(:,modes).^2,2));
    
    % Created truncated pod basis
    pod_ut  = pod_u(:,[1, modes]);
    pod_vt  = pod_v(:,[1, modes]);
    pod_vort = pod_vor(:,[1, modes]);

    if calc_coef
        
        % Calculate coefficeints for current system
        coef_problem.pod_u = pod_ut;
        coef_problem.pod_v = pod_vt;
        coef_problem.custom = custom;
        
        if custom
            fprintf('Generating coefficients for custom selected modes using %d modes\n', num_modes);
        else
            fprintf('Generating coefficients for resolved modes using %d modes\n', num_modes);
        end
        
        [l{1,1,i}, q{1,1,i}] = galerkin_coefficients(coef_problem);
        l{1,2,i} = 'Base Coeff';
        q{1,2,i} = 'Base Coeff';
        
        l{2,1,i} = galerkin_coefficients_ws(coef_problem);
        l{2,2,i} = 'Weak Coeff';
        q{2,1,i} = q{1,1,i};
        q{2,2,i} = 'Weak Coeff';
        
        % Duplicate base galerkin projection terms to all models
        l(:,:,i) = repmat(l(1:2,:,i), total_models/2, 1, 1);
        q(:,:,i) = repmat(q(1:2,:,i), total_models/2, 1, 1);
    
        %_____ Modified coefficients ______%

        % Attempt to estimate the neglected energy transfer with viscous
        fprintf('Calculating viscous dissapation terms\n');

        eddy(1,:,i) = {zeros(num_modes,1), 'Base Original'};
        eddy(2,:,i) = {zeros(num_modes,1), 'Weak Original'};

        % calculate coefficients detailed by Couplet
        if ismember('Least Squares', dissapation);
            eddy{3,1,i} = viscous_dis_couplet(modal_amp, num_modes, modes, lc{1,1}, qc{1,1}, vis);
            eddy{3,2,i} = 'Modal Base Couplet';

            eddy{4,1,i} = viscous_dis_couplet(modal_amp, num_modes, modes, lc{2,1}, qc{2,1}, vis);
            eddy{4,2,i} = 'Modal Weak Couplet';
        end

        if ismember('Averaged', dissapation);
            % Calculate coefficeints detailed by Noack
            [eddy{5,1,i}, eddy{6,1,i}] = viscous_dis(modal_amp, num_modes, lambda, l{1,1,i}, q{1,1,i}, vis);
            eddy{5,2,i} = 'Modal Base Noack';
            eddy{6,2,i} = 'Global Base Noack';

            [eddy{7,1,i}, eddy{8,1,i}] = viscous_dis(modal_amp, num_modes, lambda, l{2,1,i}, q{2,1,i}, vis);
            eddy{7,2,i} = 'Modal Weak Noack';
            eddy{8,2,i} = 'Global Weak Noack';
        end

        % Duplicate TODO name this better
        eddy(linear_models+1:total_models,1,i) = eddy(base_models+1:linear_models,1,i);
        for j = base_models+1:linear_models
            if ~isempty(eddy{j,2,i}) 
                eddy{j+linear_models-base_models,2,i} = ['NL ' eddy{j,2,i}];
            end
        end
        
        % Get correct names
        l(:,2,i) = eddy(:,2,i);
        q(:,2,i) = eddy(:,2,i);
    end
    
    %_____ Time integration _____%
    if time_int    
        
        % Load coefficient if only time integration was selected 
        if ~calc_coef 
            [l(:,:,i), q(:,:,i), eddy(:,:,i), vis] = load_coef(direct, run_num, custom);            
        end
        
        % Intial conditions and integration options
        options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
        ao = modal_amp(init,[1 modes]);
        
        % Perform final manipulation to prep integration
        [reduced_model_coeff] = integration_setup(eddy, vis, l, q, i, total_models, ...
            linear_models, num_modes);
        
        % Perform time integration
        [t, modal_amp_sim, time] = time_integration(reduced_model_coeff, eddy, vis, modal_TKE, ...
            i, t, modal_amp_sim, ao, t_scale, tspan, total_models, linear_models, options);
        disp(time);
        
        % Get correct model names
        t(:,2,i) = eddy(:,2,i);
        modal_amp_sim(:,2,i) = eddy(:,2,i);
        
        %____ Plotting functions ____%

        % Classify simulation to to empirical clusters
        if classify_sim && num_modes <= 40
            idx = (cluster_range == num_modes-1);
            [scores_km(:,1,i), scores_gm(:,1,i)] = classify_Gal(centers{idx}, ...
                gm_models{idx}, modal_amp_sim(:,:,i), direct, km_stoch{idx}, gm_stoch{idx});
        end

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
        plot_data.type          = 'Galerkin';
        plot_data.plot_type     = plot_type;
        plot_data.sample_freq   = sample_freq;
        plot_data.x             = x;
        plot_data.y             = y;
        plot_data.Mod           = false;
        plot_data.bnd_idx       = bnd_idx;
        plot_data.custom        = custom;

        all_t = t(:,1,i)';
        all_ids = eddy(:,2,i)';
        all_modal_amps =  modal_amp_sim(:,1,i)';

        close all
        
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
    
    % generate empty results structures
    results_coef = struct;
    results_int = struct;
    
    % Include coefficent terms if requested
    if calc_coef
        results_coef.modes = modes;
        results_coef.run_num = run_num;
        results_coef.num_modesG = num_modes-1;
        results_coef.sample_freq = sample_freq;
        results_coef.linear_models = linear_models;
        results_coef.total_models = total_models;

        % viscous and convective terms
        results_coef.l = squeeze(l(:,:,i));
        results_coef.q = squeeze(q(:,:,i));

        % visocity and eddy-visocity
        results_coef.eddy = squeeze(eddy(:,:,i));
        results_coef.vis = vis;
    end

    % Include integration information if requested
    if time_int
        results_int.num_modesG = num_modes-1;
        results_int.run_num = run_num;

        results_int.modal_amp_sim = squeeze(modal_amp_sim(:,:,i));
        results_int.t = squeeze(t(:,:,i));
    end
    
    % Save relavent coefficients
    if save_coef == true
        futures = save_galerkin(direct, custom, time_int, calc_coef, i, ...
            futures, results_coef, results_int);
    end
end

% Return values if requested
if nargout >= 1
    results_coef.num_modesG = num_modesG;
    res_coef = results_coef;
end
if nargout == 2
    res_int = results_int;
end

% return format
format short g
end

