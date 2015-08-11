function Mod_POD(varargin)
% MOD_POD modifies an existing Galerkin system by creating a new POD basis
% that is transformed such that it minimizing the energy decrease compared
% to the original system. will show results of new system via graphing
%
%   MOD_POD() prompt user for analysis folder for a given test run, will use
%   all defaultes detailed below
% 
%   MOD_POD(problem) Using fields provided in the structure PROBLEM sets up
%   analysis specified by PROBLEM. all unfilled fields go to defaults
%
%   problem.num_modes = 10
%   Specify the number of modes that will be calculated the basis
%   transformation
%
%   problem.basis_modes = 'double'
%   Specify the basis that you want the reduced model to be constructed 
%   from by default a basis of 2*num_modes is selected
%
%   problem.plot_type = {'amp', 'fft', 'energy'}
%   Specify which outout graphes are desired, current options are modal
%   amplitude 'amp', fourier fast transform 'fft', and 'video which produces
%   a video of the simulated flow
%
%   problem.save_mod = true
%   Save relvant values to a .mat file
%
%   problem.init = 1
%   Specify which image will constitute the initial conditions
%
%   problem.tspan = 0:0.0001:1
%   Specify time span for time integration. If tspan{1} = {"test"} then the
%   data will be integrated over the range of empirical data, if tspan{2}
%   is set to a number n it will run n times shorters
%
%   problem.line_range = 100
%   Specify a relatively range of epsilon values to search to be course
%   searched in order to generate a transformation matrix which produces an
%   energy balanced system
%
%   problem.direct = ''
%   Specify directory that will be searched for POD data, default is to
%   prompt user
%
%   problem.run_num = 'first'
%   Specify which run this Galerkn Projection should be based from, default
%   is to use the most recent
%
%   problem.models = {'GM', 'GM1'}
%   Specify which class of Galerkin models from GALKERKIN_PROJ are corrected 
%   with a basis transform. specify 'all' to use all models
%   
%   problem.submodels = {'base', 'weak'}
%   Specify which submodles will be used in 
%
%   problem.fft_window = [0 2000]
%   Specify the hertz range that the fft plot should capture
%
%   problem.score = true
%   Classify simulated results to emprical data
%
%   problem.num_cores = 'auto'
%   Set the max number of cores to be used for POD_Gen, auto will set this
%   to the number of cores in the computer.

format long g
close all
clc;

%List of fields that will be checked
fields = {  'num_modes',    'plot_type',    'save_mod', ...
            'init',         'line_range',   'direct' ,...
            'run_num',      'models',       'submodels',...   
            'fft_window',   'tspan',        'basis_modes',...
            'custom',       'score',        'num_cores', ...
            'outlier_mode'};

% Parse problem structure provided to set it up correctly
if nargin == 1
    problem = parse_inputs(fields, @setdefaults_mod, varargin{1});
else
    fprintf('Provide a single structure as input, use help Galerkin_Proj for information.\n');
    fprintf('Using Defaults\n\n');
    problem = parse_inputs(fields, @setdefaults_mod);
end

num_modes   = problem.num_modes;
basis_modes = problem.basis_modes;
plot_type   = problem.plot_type;
save_mod    = problem.save_mod;
init        = problem.init;
line_range  = problem.line_range;
direct      = problem.direct;
run_num     = problem.run_num;
models      = problem.models;
submodels   = problem.submodels;
fft_window  = problem.fft_window;
tspan       = problem.tspan;
custom      = problem.custom;
score       = problem.score;
num_cores   = problem.num_cores;
outlier_mode    = problem.outlier_mode;


clear problem

% Setup MATLAB to a max number of cores
setup_cores(num_cores);

% Handle File IO
if strcmp(direct, '');
    [direct_POD, direct] = prompt_folder('POD', run_num);
    [direct_Gal, direct] = prompt_folder('Galerkin', run_num, direct, basis_modes, custom);
else
    [direct_POD, direct] = prompt_folder('POD', run_num, direct);
    [direct_Gal, direct] = prompt_folder('Galerkin', run_num, direct, basis_modes, custom);
end


% Make sure folders are up to date and load collected data
update_folders(direct);
fix_gmm(direct);

% Load POD variables
vars = load(direct_POD, 'results_pod');

% Create more readable names
pod_U       = vars.results_pod.pod_U;       % pod modes
mean_U      = vars.results_pod.mean_U;      % mean flow
u_scale     = vars.results_pod.u_scale;     % velocity scaling
l_scale     = vars.results_pod.l_scale;     % length scaling
lambda_basis= vars.results_pod.lambda;      % eigenvalues of modes
non_dim     = vars.results_pod.non_dim;     % was the result non_dimensionalized
modal_amp   = vars.results_pod.modal_amp;   % empirical modal amp
dimensions  = vars.results_pod.dimensions;  % dimensions of mesh
run_num     = vars.results_pod.run_num;     % POD run numbers
X           = vars.results_pod.X;           % x coordinate
bnd_idx     = vars.results_pod.bnd_idx;     % boundary matrix
exp_sampling_rate = vars.results_pod.exp_sampling_rate; % empirical sampling

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

% Create non-dimensionalized timescale
[t_scale, tspan] = calc_t_scale(u_scale, l_scale, non_dim, tspan);     

% If requested load cluster data
if score
    
    % Load Cluster variables
    vars = load(direct_POD, 'results_clust');
    
    km = vars.results_clust.km;     % k-means cluster data
    gmm = vars.results_clust.gmm;     % gaussian mixture data
    num_clusters = vars.results_clust.num_clusters;

    
    frob_km     = cell(length(num_modes),1);
    frob_gmm    = cell(length(num_modes),1);
    prob_km     = cell(length(num_modes),1);    
    prob_gmm    = cell(length(num_modes),1);   
    completed   = cell(length(num_modes),1);
end

vars = load(direct_Gal, 'results_coef');

% Get results from Galerkin Projection
system      = vars.results_coef.system;
modes       = vars.results_coef.modes;

if ischar(models) && strcmp(models, 'all')
    models = flow_comps(system.eddy);
    m = models;
else
    m = models;
end

if ischar(submodels) && strcmp(submodels, 'all')
    submodels = flow_comps(system.eddy.(m{1}));
    s = submodels;
else
    s = submodels;
end

% Free memory
clear vars

[~, u] = flow_comps(X, pod_U);
dims = flow_ncomps(X);

% Main Loop
for i = 1:length(models)
    for j = 1:length(submodels)
                
        % Set variables for model
        [C, L, Q, lambda, modal_amp] = term2order(system, ...
                    lambda_basis, modal_amp, m{i}, s{j}, modes);
        
        % Create temporary pod basis from requested modes
        for k = 1:dims
            pod_Ut.(u{k}) = pod_U.(u{k})(:,modes);
        end
        modal_ampt = modal_amp(:, modes);
        
        
        % Creat line search problem
        line_problem.C = C;
        line_problem.L = L;
        line_problem.Q = Q;
        line_problem.t_scale = t_scale;
        line_problem.tspan = tspan;
        line_problem.init = init;
        line_problem.basis_modes = basis_modes;
        line_problem.num_modes = num_modes;
        line_problem.lambda = lambda;
        line_problem.modal_amp = modal_ampt;
        line_problem.line_range = line_range;

        
%         % Brute force line search
        [epsilon_low, epsilon_high, ~, flip] = line_search(line_problem);
        
        if flip == true
            epsilon_range = [epsilon_low, epsilon_high];
            options = optimset('PlotFcns', {@optimplotx, @optimplotfval}, 'Display', 'iter', ...
                'FunValCheck', 'on');
            
            [epsilon_final, ~, ~, OUTPUT] = fzero(@(epsilon) optimal_rotation...
                (epsilon, C, L, Q, num_modes, basis_modes, lambda, modal_ampt, t_scale, tspan, init, 64000), epsilon_range, options);
            disp(OUTPUT);
        else
            disp('no sign flip detected');
            options = optimset('PlotFcns', {@optimplotx, @optimplotfval}, 'Display', 'iter', 'TolFun', 1e-10, 'TolX', 1e-6);
            
            [epsilon_final, ~, ~, OUTPUT] = fminbnd(@(epsilon) abs(optimal_rotation...
                (epsilon, C, L, Q, num_modes, basis_modes, lambda, modal_ampt, t_scale, tspan, init, 64000)), epsilon_low, epsilon_high, options);
            disp(OUTPUT);
        end
        
        close all;
        
        s_til{j} = ['til_' s{j}];
        
        % Final calculation of transformation matrix and new constant linear and
        % quadratic terms
        [~, X_til, C_til, L_til, Q_til, modal_amp_til, t] = ...
            optimal_rotation(epsilon_final, C, L, Q, num_modes, basis_modes, lambda, modal_ampt, t_scale, tspan, init, 64000);
        
        [pod_U_til, modal_amp_raw_til] = ...
            basis_transform(pod_Ut, modal_ampt, num_modes, X_til);
        
        integration.t.(m{i}).(s_til{j}) = t;
        integration.modal_amp.(m{i}).(s_til{j}) = modal_amp_til;
        
        % Prepare data to be saved
        system.X_til.(m{i}).(s_til{j}) = X_til;
        system.C_til.(m{i}).(s_til{j}) = C_til;
        system.L_til.(m{i}).(s_til{j}) = L_til;
        system.Q_til.(m{i}).(s_til{j}) = Q_til;
        system.X_til.(m{i}).(s_til{j}) = X_til;
        system.pod_U_til.(m{i}).(s_til{j}) = pod_U_til;
        system.modal_amp_til.(m{i}).(s_til{j}) = modal_amp_raw_til;
        system.epsilon.(m{i}).(s_til{j}) = epsilon_final;
        
         % Classify simulation to to empirical clusters
        if score 
            
            % Ready score info Structure
            score_info.km = km;
            score_info.gmm = gmm;
            score_info.integration = integration;
            score_info.tspan = tspan;
            score_info.modal_amp = modal_amp_raw_til;
            score_info.num_clusters = num_clusters;
            score_info.int_time = 3600;
            score_info.num_cores = num_cores;
            score_info.modes = 1:num_modes;
            score_info.custom = true;
            score_info.outlier_mode = outlier_mode;
            score_info.direct = direct;
            score_info.multiplier = multiplier;
            score_info.MOD = true;
            score_info.modal_amp_til = system.modal_amp_til;
            
            % score results
            [frob_km{i}, frob_gmm{i}, prob_km{i}, prob_gmm{i}, completed{i}] = ...
                score_model(score_info);
        end
        
        % Prepare data
        plot_data.num_modes     = num_modes;
        plot_data.direct        = direct;
        plot_data.pod_Ut        = pod_U_til;
        plot_data.mean_U        = mean_U;
        plot_data.dimensions    = dimensions;
        plot_data.fft_window    = fft_window;
        plot_data.u_scale       = u_scale;
        plot_data.l_scale       = l_scale;
        plot_data.type          = 'MOD';
        plot_data.Mod           = true;
        plot_data.plot_type     = plot_type;
        plot_data.sample_freq   = sample_freq;
        plot_data.id            = strrep([m{i} ' ' s_til{j}], '_', ' ');
        plot_data.modal_amp     = modal_amp_til;
        plot_data.t             = t;
        plot_data.X             = X;
        plot_data.custom        = false;
        plot_data.bnd_idx       = bnd_idx;
        
        % Generate plots
        produce_plots(plot_data);
       
    end
end

results_mod_coef.name           = 'results_mod_coef';
results_mod_int.name            = 'results_mod_int';
results_mod_scores.name         = 'results_mod_scores';

results_mod_coef.run_num        = run_num;
results_mod_coef.modes          = 1:num_modes;
results_mod_coef.system         = system;

results_mod_int.integration     = integration;
results_mod_int.run_num         = run_num;

if score
    results_mod_scores.frob_km      = frob_km;
    results_mod_scores.frob_gmm     = frob_gmm;
    results_mod_scores.prob_km      = prob_km;
    results_mod_scores.prob_gmm     = prob_gmm;
    results_mod_scores.completed    = completed;
    results_mod_scores.run_num      = run_num;
end

if save_mod == true
    save_results(num_modes, direct, 'Mod Galerkin Coeff', custom, results_mod_coef, ...
        results_mod_int, results_mod_scores);
end

end



