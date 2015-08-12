function [res_coef, res_int, res_scores] = Galerkin_Proj(varargin)
% GALERKIN_PROJ perform Galerkin projection of the POD system onto Navier
% Stokes. This requires POD_GEN to be run, can produce and save graphs
% see below return results into RES if output is requested
%
%   GALERKIN_PROJ() prompt user for analysis folders for a given test run,
%   will use all defaults detailed below
% 
%   GALERKIN_PROJ(problem) Using fields provided in the structure PROBLEM
%   sets up analysis specified by PROBLEM. all unfilled fields go to defaults
%
%   problem.num_modesG = 10
%   Specify the number of modes that will calculated in Galerkin projection
%
%   problem.plot_type = {'amp', 'fft', 'energy'}
%   Specify which outout graphes are desired, current options are modal
%   amplitude 'amp', fourier fast transform 'fft', system energy 'energy', 
%   'video' which produces a video of the simulated flow using quivers, 
%   'video stream' will produce a video of simulated flow using streamlines,
%   caution very slow
% 
%   problem.save_coef = true
%   Save relvant values to at .mat file
%
%   problem.override_coef = false
%   Overwrite previous coefficients for the same run number, if the run
%   numbers are different new coefficients will be calcfinteulated
%
%   problem.tspan = {'test'}
%   Specify time span for time integration. If tspan{1} = {"test"} then the
%   data will be integrated over the range of empirical data, if tspan{2}
%   is set to a number n it will run n times shorters
%
%   problem.int_time = 3600
%   Specifiy the maximum integration time in seconds.
% 
%   problem.init = 1
%   Specify which image will constitute the initial conditions
%
%   problem.direct = ''
%   Specify directory that will be searched for POD data, default is to
%   prompt user
%
%   problem.Re0_gen = @Re0_gen_shear
%   Fuction handle used to produced Reynolds number for the run, takes one
%   input of working project directory
%
%   problem.fft_window = [0 2000]
%   Specify the hertz range that the fft plot should capture
%
%   problem.run_num = 'first'
%   Specify which run this Galerkn Projection should be based from, default
%   is to use the most recent
%
%   problem.dissapation = {'Least Squares', 'Averaged'}
%   Specify which viscious dissapation method(s) to use
%
%   problem.time_int = true
%   If false no time integration or plotting will be performed, used to
%   generated the Galerkin systems for MOD_POD
%
%   problem.calc_coef = true
%   If false attempt to locate previous runs files to time integration
%
%   problem.use_chunks = false
%   Set this to true if you are running out of memeory well write values of q
%   to disk 
%
%   problem.score = true
%   Score produces models using Surrogate Markov Models
% 
%   problem.odesolver = @ode113
%   Select the ode solver to be used for time integration.
%
%   problem.num_cores = 'auto'
%   Set the max number of cores to be used for POD_Gen, auto will set this
%   to the number of cores in the computer.

% Set format, clear figures, and set up correct directory
format long g
close all
clc;

%List of fields that will be checked
fields = {  'num_modesG',   'plot_type',    'save_coef', ...
            'override_coef','tspan',        'init', ...
            'direct' ,      'Re0_gen',      'fft_window', ...
            'run_num',      'dissapation',  'time_int', ...
            'use_chunks',   'calc_coef',    'score', ...
            'odesolver',    'int_time',     'num_cores', ...
            'outlier_mode'};
        
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
int_time        = problem.int_time;
tspan           = problem.tspan;
init            = problem.init;
direct          = problem.direct;
Re0_gen         = problem.Re0_gen;     
fft_window      = problem.fft_window;
dissapation     = problem.dissapation;
calc_coef       = problem.calc_coef;
time_int        = problem.time_int;
score           = problem.score;
use_chunks      = problem.use_chunks;
odesolver       = problem.odesolver;
num_cores       = problem.num_cores;
outlier_mode    = problem.outlier_mode;


clear problem

% Setup MATLAB to a max number of cores
setup_cores(num_cores);

fprintf('\nLoading POD variables\n\n');

% Prompt User for folder if directory is not provided
if strcmp(direct, '');
    [direct_POD, direct] = prompt_folder('POD', run_num);
else
    [direct_POD, direct] = prompt_folder('POD', run_num, direct);
end

% Check folders are up to most recent format
update_folders(direct);
fix_gmm(direct);

% Load POD variables
vars = load(direct_POD);

% Create more readable names
X           = vars.results_pod.X;           % mesh coordinates
u_scale     = vars.results_pod.u_scale;     % velocity scaling
l_scale     = vars.results_pod.l_scale;     % length scaling
pod_U       = vars.results_pod.pod_U;       % pod modes
lambda      = vars.results_pod.lambda;      % eigenvalues of modes
dimensions  = vars.results_pod.dimensions;  % dimensions of mesh
bnd_idx     = vars.results_pod.bnd_idx;     % location of boundaries
non_dim     = vars.results_pod.non_dim;     % was the result non_dimensionalized
run_num     = vars.results_pod.run_num;     % POD run numbers
cutoff      = vars.results_pod.cutoff;      % number of modes at cutoff
modal_amp   = vars.results_pod.modal_amp;   % modal amplitude  from raw data
exp_sampling_rate = vars.results_pod.exp_sampling_rate;

if calc_coef
    bnd_X       = vars.results_pod.bnd_X;     % location of flow boundaries normal to x
    volume      = vars.results_pod.volume;    % mesh area size
end

if time_int && score
%     if isfield(vars, 'results_clust');
%         km = vars.results_clust.km; % k-means cluster information
%         gmm = vars.results_clust.gmm; % gaussian mixture model cluster information
%         num_clusters = vars.results_clust.num_clusters; % number of clustered used
%     else
        % Fill with stub if empty
        km = struct;
        gmm = struct;
        num_clusters = 10;
%     end
end

% Free memory
clear vars 


% Initialize coefficeints
system    = cell(length(num_modesG),1);

% Determine Reynolds number 
Re0 = Re0_gen(direct);  

% Determine kinematic visocity
for i = 1:length(num_modesG)
    system{i}.vis = calc_vis(u_scale, l_scale, Re0, non_dim);
end

% Determine sampling frequency from provided tspan
if isnumeric(tspan)
    sample_freq = 1/(tspan(2) - tspan(1));
    multiplier = 1;
    fprintf('Simulated Sampling Frequency %6.2f Hz\n\n', sample_freq);
end
if iscell(tspan) && strcmp(tspan{1}, 'test')
    [tspan, sample_freq, multiplier] = calc_tspan(tspan, exp_sampling_rate, modal_amp);
    fprintf('Simulated Sampling Frequency %6.2f Hz\n\n', sample_freq);
end

% Create modify time scale if non-dimenionalized
[t_scale, tspan] = calc_t_scale(u_scale, l_scale, non_dim, tspan);

% Preallocate for future save function
futures = parallel.FevalOnAllFuture;

%% Setup variables to Coefficient Calculations

if calc_coef 
    
    % Ready Coef Problem Structure
    coef_problem.X              = X;
    coef_problem.use_chunks     = use_chunks;
    coef_problem.pod_U          = pod_U;
    coef_problem.dimensions     = dimensions;
    coef_problem.volume         = volume;
    coef_problem.bnd_idx        = bnd_idx;
    coef_problem.bnd_X          = bnd_X;
    coef_problem.run_num        = run_num;
    coef_problem.override_coef  = override_coef;
    coef_problem.direct         = direct;
    coef_problem.custom         = false;
    coef_problem.num_cores      = num_cores;
    
    % Free memory
    clear volume bnd_X 
    
    % Generate values up to cutoff to be used for Couplet viscous dissipation
    if ismember('Least Squares', dissapation);
        fprintf('Generating coefficients for unresolved modes using %d modes\n', cutoff);
        
        % Calculate for traditional Galerkin Projection
        [cutoff_coef.l.base, cutoff_coef.q.base] = galerkin_coefficients(coef_problem);
        
        % Calculate for weak formulation Galerkin Projection
        cutoff_coef.l.weak  = galerkin_coefficients_ws(coef_problem);
        cutoff_coef.q.weak = cutoff_coef.q.base;
    end 
end

[~, u] = flow_comps_ip(X, pod_U);
dims = flow_dims(X);

% Initialize if time integration is being performed
if time_int
    integration = cell(length(num_modesG), 1);
end

% Initialize scores if clustering
if score
    frob_km     = cell(length(num_modesG), 1);
    frob_gmm     = cell(length(num_modesG), 1);
    like_km     = cell(length(num_modesG), 1);    
    like_gmm     = cell(length(num_modesG), 1);
    completed   = cell(length(num_modesG), 1);
end

% Perform calculation on each set of modes requested
for i = 1:length(num_modesG)
    fprintf('\n');
    
    integration{i}.t = struct;
    integration{i}.modal_amp = struct;
    
    % Determine if custom modes have been used, and setup appropriately 
    if iscell(num_modesG)
        custom = true;
        num_modes = length(num_modesG{i}) + 1;
        modes = num_modesG{i} + 1;
    else
        custom = false;
        num_modes = num_modesG(i)+1;
        modes = 2:num_modes;
    end
    
    % Calculate empirical average TKE for selected modes
    modal_TKE = mean(sum(1/2*modal_amp(:,modes).^2,2));
    
    % Create temporary pod basis from requested modes
    for j = 1:dims
        pod_Ut.(u{j}) = pod_U.(u{j})(:,[1, modes]);
    end

    if calc_coef
        
        % Calculate coefficeints for current system
        coef_problem.pod_U = pod_Ut;
        coef_problem.custom = custom;
        
        if custom
            fprintf('Generating coefficients for custom selected modes using %d modes\n', num_modes);
        else
            fprintf('Generating coefficients for resolved modes using %d modes\n', num_modes);
        end
        
        
        [system{i}.l.base, system{i}.q.base] = galerkin_coefficients(coef_problem);
        
        system{i}.l.weak = galerkin_coefficients_ws(coef_problem);
        system{i}.q.weak = system{i}.q.base; 
        
        s = flow_comps(system{i}.l);
        s_comps = flow_ncomps(system{i}.l);
        
        %_____ Modified coefficients ______%

        % Attempt to estimate the neglected energy transfer with viscous
        fprintf('Calculating viscous dissapation terms\n');

        for j = 1:s_comps
            system{i}.eddy.GM.(s{j}) = zeros(num_modes,1);
        end

        % Averaged energy balance detailed by Noack
        if ismember('Averaged', dissapation);
            for j = 1:s_comps
            [system{i}.eddy.GM2.(s{j}), system{i}.eddy.GM1.(s{j})] = viscous_dis(modal_amp, ...
                num_modes, lambda, system{i}.l.(s{j}), system{i}.q.(s{j}), system{i}.vis);
            end
        end
        
        % Energy balance via least squares solution detailed by Couplet
        if ismember('Least Squares', dissapation);
            for j = 1:s_comps
            system{i}.eddy.GM3.(s{j}) = viscous_dis_couplet(modal_amp, num_modes, modes, ...
                            cutoff_coef.l.(s{j}), cutoff_coef.q.(s{j}), system{i}.vis);
            end
        end
    end
    
    %_____ Time integration _____%
    if time_int    
        
        % Load coefficient if only time integration was selected 
        if ~calc_coef 
            system{i} = load_coef(direct, run_num, num_modes, custom);            
        end
        
        % Intial conditions and integration options
        options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
        ao = modal_amp(init,[1 modes]);
        
        % Perform final manipulation to prep integration
        system = integration_setup(system, i, num_modes);
        
        % Perform time integration
        [integration, time] = time_integration(odesolver, system, integration, modal_TKE, ...
                                i, ao, t_scale, tspan, int_time, options);
        disp(time);
        
        %____ Plotting functions ____%

        % Classify simulation to to empirical clusters
        if score
            
            % Ready score info Structure
            score_info.km = km;
            score_info.gmm = gmm;
            score_info.integration = integration{i};
            score_info.tspan = tspan;
            score_info.modal_amp = modal_amp;
            score_info.num_clusters = num_clusters;
            score_info.int_time = int_time;
            score_info.num_cores = num_cores;
            score_info.modes = modes;
            score_info.custom = custom;
            score_info.outlier_mode = outlier_mode;
            score_info.direct = direct;
            score_info.multiplier = multiplier;
            score_info.MOD = false;
            
            % score results
            [frob_km{i}, frob_gmm{i}, like_km{i}, like_gmm{i}, completed{i}] = ...
                score_model(score_info);
        end

        % Prepare data
        plot_data.num_modes     = num_modes-1;
        plot_data.direct        = direct;
        plot_data.pod_Ut        = pod_Ut;
        plot_data.dimensions    = dimensions;
        plot_data.fft_window    = fft_window;
        plot_data.u_scale       = u_scale;
        plot_data.l_scale       = l_scale;
        plot_data.type          = 'Galerkin';
        plot_data.plot_type     = plot_type;
        plot_data.sample_freq   = sample_freq;
        plot_data.X             = X;
        plot_data.Mod           = false;
        plot_data.bnd_idx       = bnd_idx;
        plot_data.custom        = custom;

        close all
        
        m = flow_comps(integration{i}.t);
        m_comps = flow_ncomps(integration{i}.t);
        
        % Cycle through and plot all requested figures
        for j = 1:m_comps
            s = flow_comps(integration{i}.t.(m{j}));
            s_comps = flow_ncomps(integration{i}.t.(m{j}));
            for k = 1:s_comps
                if isempty(integration{i}.t.(m{j}).(s{k})) || ...
                        isfield(integration{i}.t.(m{j}).(s{k}), 'error')
                    continue;
                end
                plot_data.t = integration{i}.t.(m{j}).(s{k});
                plot_data.id = strrep([m{j} ' ' s{k}], '_', ' ');
                plot_data.modal_amp = integration{i}.modal_amp.(m{j}).(s{k});
                produce_plots(plot_data);
            end
        end
    end
        
    fprintf('Saving Galerkin Variables\n');
    
    % generate empty results structures
    results_coef = struct;
    results_int = struct;
    results_scores = struct;
    
    % Include coefficent terms if requested
    if calc_coef
        results_coef.name = 'results_coef';
        results_coef.modes = modes;
        results_coef.run_num = run_num;
        results_coef.num_modesG = num_modes-1;
        results_coef.sample_freq = sample_freq;

        % System coefficients and properties
        results_coef.system = system{i};
    end

    % Include integration information if requested
    if time_int
        results_int.name = 'results_int';
        results_int.num_modesG = num_modes-1;
        results_int.run_num = run_num;
        results_int.tspan = tspan;

        results_int.integration = integration{i};
    end
    
    % Include scores if requested 
    if score
        results_scores.name = 'results_scores';
        results_scores.frob_km = frob_km{i};
        results_scores.frob_gmm = frob_gmm{i};
        results_scores.like_km = like_km{i};
        results_scores.like_gmm = like_gmm{i};
        results_scores.completed = completed{i};
    end
    
    % Save relavent coefficients
    if save_coef == true
        futures = save_galerkin(direct, custom, time_int, calc_coef, score, i, ...
            futures, num_modes-1, results_coef, results_int, results_scores);
    end
end

% Return values if requested
if nargout >= 1
    results_coef.num_modesG = num_modesG;
    res_coef = results_coef;
end
if nargout >= 2
    res_int = results_int;
end
if nargout == 3
    res_scores = results_scores;
end


% return format
format short g
end

