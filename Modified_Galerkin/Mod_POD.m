function Mod_POD(varargin)
% MOD_POD modifies an existing Galerkin system by creating a new POD basis
% that is transformed such that it minimizing the energy decrease compared
% to the original system. will show results of new system via graphing
%
% MOD_POD() prompt user for analysis folder for a given test run, will use
% all defaultes detailed below
% 
% MOD_POD(problem) Using fields provided in the structure PROBLEM sets up
% analysis specified by PROBLEM. all unfilled fields go to defaults
%
% problem.RD_nm = 10
% Specify the number of modes that will be calculated the basis
% transformation
%
% problem.OG_nm = 'double'
% Specify the basis that you want the reduced model to be constructed from
% by default a basis of 2*RD_nm is selected
%
% problem.plot_type = {'amp', 'fft', 'energy'}
% Specify which outout graphes are desired, current options are modal
% amplitude 'amp', fourier fast transform 'fft', and 'video which produces
% a video of the simulated flow
%
% problem.save_mod = true
% Save relvant values to a .mat file
%
% problem.init = 1
% Specify which image will constitute the initial conditions
%
% problem.tspan = 0:0.0001:1
% Specify for how long integration will be performed will be sampled at
% sampling frequency of base data
%
% problem.line_range = 100
% Specify a relatively range of epsilon values to search to be course
% searched in order to generate a transformation matrix which produces an
% energy balanced system
%
% problem.direct = ''
% Specify directory that will be searched for POD data, default is to
% prompt user
%
% problem.run_num = 'first'
% Specify which run this Galerkn Projection should be based from, default
% is to use the most recent
%
% problem.models = 1:8
% Specify which Galerkin Models from GALKERKIN_PROJ are corrected with a
% basis transform. 
%
% problem.fft_window = [0 2000]
% Specify the hertz range that the fft plot should capture

format long g
close all
clc;

%List of fields that will be checked
fields = {  'RD_nm',        'plot_type',    'save_mod', ...
            'init',         'line_range',   'direct' ,...
            'run_num',      'models',       'fft_window', ...
            'tspan',        'OG_nm'         'custom'};

% Parse problem structure provided to set it up correctly
if nargin == 1
    problem = parse_inputs(fields, @setdefaults_mod, varargin{1});
else
    fprintf('Provide a single structure as input, use help Galerkin_Proj for information.\n');
    fprintf('Using Defaults\n\n');
    problem = parse_inputs(fields, @setdefaults_mod);
end

RD_nm       = problem.RD_nm;
OG_nm       = problem.OG_nm;
plot_type   = problem.plot_type;
save_mod    = problem.save_mod;
init        = problem.init;
line_range  = problem.line_range;
direct      = problem.direct;
run_num     = problem.run_num;
models      = problem.models;
fft_window  = problem.fft_window;
tspan       = problem.tspan;
custom      = problem.custom;

% Check status of parrallel pool
if isempty(gcp('nocreate'));
    parpool('local');
end

% Handle File IO
if strcmp(direct, '');
    [direct_POD, direct] = prompt_folder('POD', run_num);
    [direct_Gal, direct] = prompt_folder('Galerkin', run_num, direct, OG_nm, custom);
else
    [direct_POD, direct] = prompt_folder('POD', run_num, direct);
    [direct_Gal, direct] = prompt_folder('Galerkin', run_num, direct, OG_nm, custom);
end

% Make sure folders are up to date and load collected data
update_folders(direct);

% Load POD variables
vars = load(direct_POD, 'results_pod');

% Create more readable names
pod_u       = vars.results_pod.pod_u;       % streamwise pod modes
pod_v       = vars.results_pod.pod_v;       % spanwise pod modes
mean_u      = vars.results_pod.mean_u;
mean_v      = vars.results_pod.mean_v;
u_scale     = vars.results_pod.u_scale;     % velocity scaling
l_scale     = vars.results_pod.l_scale;     % length scaling
lambda_OG   = vars.results_pod.lambda;      % eigenvalues of modes
modal_amp   = vars.results_pod.modal_amp;   % modal amplitude of each image in pod basis
dimensions  = vars.results_pod.dimensions;  % dimensions of mesh
run_num     = vars.results_pod.run_num;     % POD run numbers
x           = vars.results_pod.x;           % x coordinate
y           = vars.results_pod.y;           % y coordinate
bnd_idx     = vars.results_pod.bnd_idx;

% Determine sampling frequency from provided tspan
if length(tspan) > 2
    sample_freq = 1/(tspan(2) - tspan(1));
    fprintf('Detected Sampling Frequency %6.2f Hz\n\n', sample_freq);
else
    error('must provide tspan with a range');
end

% Create non-dimensionalized timescale
t_scale = u_scale/l_scale;  
tspan = tspan*t_scale;      

vars = load(direct_POD, 'results_clust');
    
centers         = vars.results_clust.centers;     % k-means cluster centers
km_stoch        = vars.results_clust.km_stoch;    % k-mean stochastic matrix
gm_stoch        = vars.results_clust.gm_stoch;    % gaussian mixture stochastic matrix
gm_models       = vars.results_clust.gm_models;   % gaussian mixture models
cluster_range   = vars.results_clust.cluster_range;    % number of variables in cluster

vars = load(direct_Gal, 'results_coef');

% Get results from Galerkin Projection
q_total     = vars.results_coef.q;
l_total     = vars.results_coef.l;
eddy_total  = vars.results_coef.eddy;
vis         = vars.results_coef.vis;
linear_models = vars.results_coef.linear_models;

% Free memory
clear vars

% Remove any passed models that are greater than the linear models
models(models > linear_models) = [];

% Main Loop
for i = 1:size(models, 2)
    
    % Calculate total visocity
    total_vis = vis + eddy_total{models(i),1};
    
    % Set variables for model
    [C, L, Q, lambda, modal_amp] = term2order(l_total{models(i),1}, q_total{models(i),1}, ...
                                   total_vis, lambda_OG, modal_amp, modes);
    
    % Truncate POD
    pod_ut = pod_u(:,modes);
    pod_vt = pod_v(:,modes);
    modal_ampt = modal_amp(:, modes);
    

    % Creat line search problem
    line_problem.C = C;
    line_problem.L = L;
    line_problem.Q = Q;
    line_problem.t_scale = t_scale;
    line_problem.tspan = tspan;
    line_problem.init = init;
    line_problem.OG_nm = OG_nm;
    line_problem.RD_nm = RD_nm;
    line_problem.lambda = lambda;
    line_problem.modal_amp = modal_ampt;
    line_problem.line_range = line_range;

    % Brute force line search
    [epsilon_low, epsilon_high, ~, flip] = line_search(line_problem);

    if flip == true
        epsilon_range = [epsilon_low, epsilon_high];
        options = optimset('PlotFcns', {@optimplotx, @optimplotfval}, 'Display', 'iter', ...
        'FunValCheck', 'on');

        [epsilon_final, ~, ~, OUTPUT] = fzero(@(epsilon) optimal_rotation...
            (epsilon, C, L, Q, OG_nm, RD_nm, lambda, modal_ampt, t_scale, tspan, init, 64000), epsilon_range, options);
        disp(OUTPUT);
    else
        disp('no sign flip detected');
        options = optimset('PlotFcns', {@optimplotx, @optimplotfval}, 'Display', 'iter', 'TolFun', 1e-10, 'TolX', 1e-6);

        [epsilon_final, ~, ~, OUTPUT] = fminbnd(@(epsilon) abs(optimal_rotation...
            (epsilon, C, L, Q, OG_nm, RD_nm, lambda, modal_ampt, t_scale, tspan, init, 64000)), epsilon_low, epsilon_high, options);
        disp(OUTPUT);
    end
    
    close all;
    

    % Final calculation of transformation matrix and new constant linear and
    % quadratic terms
    [~, X, C_til, L_til, Q_til, modal_amp_til, t] = ...
        optimal_rotation(epsilon_final, C, L, Q, OG_nm, RD_nm, lambda, modal_ampt, t_scale, tspan, init, 64000);

    [pod_u_til, pod_v_til, modal_amp_raw_til] = ...
        basis_transform(pod_ut, pod_vt, modal_ampt, RD_nm, X);
    
%     idx = (cluster_range == RD_nm);
%     [scores_km, scores_gm] = classify_Gal(centers{idx}, ...
%         gm_models{idx}, modal_amp_til, direct, km_stoch{idx}, gm_stoch{idx});

    % Prepare data
    plot_data.num_modes     = RD_nm;
    plot_data.direct        = direct;
    plot_data.init          = init;
    plot_data.pod_ut        = pod_u_til;
    plot_data.pod_vt        = pod_v_til;
    plot_data.mean_u        = mean_u;
    plot_data.mean_v        = mean_v;
    plot_data.dimensions    = dimensions;
    plot_data.fft_window    = fft_window;
    plot_data.u_scale       = u_scale;
    plot_data.l_scale       = l_scale;
    plot_data.plot_type     = plot_type;
    plot_data.sample_freq   = sample_freq;
    plot_data.type          = 'MOD';
    plot_data.id            = eddy_total{models(i),2};
    plot_data.modal_amp     = modal_amp_til;
    plot_data.t             = t;
    plot_data.x             = x;
    plot_data.y             = y;
    plot_data.Mod           = true;
    plot_data.custom        = false;
    plot_data.bnd_idx       = bnd_idx;

    % Generate plots
    produce_plots(plot_data);

    results_mod_coef.X               = X;
    results_mod_coef.C_til           = C_til;
    results_mod_coef.L_til           = L_til;
    results_mod_coef.Q_til           = Q_til;
    results_mod_coef.pod_u_til       = pod_u_til;
    results_mod_coef.pod_v_til       = pod_v_til;
    results_mod_coef.modal_amp_til   = modal_amp_raw_til;
    results_mod_coef.epsilon_final   = epsilon_final;
    
    results_mod_int.t = t;
    results_mod_int.modal_amp_til = modal_amp_til;
%     results.scores_km       = scores_km;
%     results.scores_gm       = scores_gm;
    
    % Save important coefficients
    if save_mod == true
        save([direct '\Mod Galerkin Coeff\Coeff_' num2str(run_num) '_m' num2str(RD_nm) '_' eddy_total{models(i),2} '.mat'], ...
            'results_mod_coef', 'results_mod_int');
    end
end

end



