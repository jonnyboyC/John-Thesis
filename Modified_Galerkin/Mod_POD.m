function Mod_POD(varargin)
% TODO fill out help 

format long g
close all
clc;

%List of fields that will be checked
fields = {  'RD_nm',        'plot_type',    'save_mod', ...
            'init',         'line_range',   'direct' ,...
            'run_num',      'models',       'fft_window'};

% Parse problem structure provided to set it up correctly
if nargin == 1
    problem = parse_inputs(fields, @setdefaults_mod, varargin{1});
else
    fprintf('Provide a single structure as input, use help Galerkin_Proj for information.\n');
    fprintf('Using Defaults\n\n');
    problem = parse_inputs(fields, @setdefaults_mod);
end

RD_nm       = problem.RD_nm;
plot_type   = problem.plot_type;
save_mod    = problem.save_mod;
init        = problem.init;
line_range  = problem.line_range;
direct      = problem.direct;
run_num     = problem.run_num;
models      = problem.models;
fft_window  = problem.fft_window;

% Check status of parrallel pool
if isempty(gcp('nocreate'));
    parpool('local');
end

% Handle File IO
if strcmp(direct, '');
    [direct_POD, direct] = prompt_folder('POD', run_num);
    [direct_Gal, direct] = prompt_folder('Galerkin', run_num, direct, 2*RD_nm);
else
    [direct_POD, direct] = prompt_folder('POD', run_num, direct);
    [direct_Gal, direct] = prompt_folder('Galerkin', run_num, direct, 2*RD_nm);
end

% Make sure folders are up to date and load collected data
update_folders(direct);

% TODO load variables based on type provided
vars = load(direct_POD, 'results');

% Create more readable names
pod_u       = vars.results.pod_u;       % streamwise pod modes
pod_v       = vars.results.pod_v;       % spanwise pod modes
pod_vor     = vars.results.pod_vor;     % vorticity pod modes
u_scale     = vars.results.u_scale;     % velocity scaling
l_scale     = vars.results.l_scale;     % length scaling
lambda_OG   = vars.results.lambda;      % eigenvalues of modes
modal_amp   = vars.results.modal_amp;   % modal amplitude of each image in pod basis
dimensions  = vars.results.dimensions;  % dimensions of mesh
run_num     = vars.results.run_num;     % POD run numbers
x           = vars.results.x;
y           = vars.results.y;
bnd_idx     = vars.results.bnd_idx;

clear vars

vars = load(direct_Gal, 'results');

% Get results from Galerkin Projection
q_total     = vars.results.q;
l_total     = vars.results.l;
t_total     = vars.results.t;
eddy_total  = vars.results.eddy;
vis         = vars.results.vis;
modal_amp_sim_total = vars.results.modal_amp_sim;
linear_models = vars.results.linear_models;
num_modesG = vars.results.num_modesG;

clear vars

% Remove any passed models that are greater than the linear models
models(models > linear_models) = [];

% Main Loop
for i = 1:size(models, 2)
    
    % Calculate total visocity
    total_vis = vis + eddy_total{models(i),1};
    
    % Set variables for model
    [C, L, Q, lambda, OG_nm] = term2order(l_total{models(i),1}, q_total{models(i),1}, ...
                                   total_vis, lambda_OG, num_modesG);
    t = t_total{models(i),1};
    modal_amp_sim = modal_amp_sim_total{models(i),1};
    

    % Truncate POD
    pod_ut = pod_u(:,1:OG_nm);
    pod_vt = pod_v(:,1:OG_nm);
    pod_vort = pod_vor(:,1:OG_nm);
    modal_ampt = modal_amp(:, 1:OG_nm);
    

    line_problem.C = C;
    line_problem.L = L;
    line_problem.Q = Q;
    line_problem.t = t;
    line_problem.init = init;
    line_problem.OG_nm = OG_nm;
    line_problem.RD_nm = RD_nm;
    line_problem.lambda = lambda;
    line_problem.modal_amp = modal_ampt;
    line_problem.line_range = line_range;

    [epsilon, transfer, flip_idx] = line_search(line_problem);

    if flip_idx ~= 0
        epsilon_range = [epsilon(flip_idx) epsilon(flip_idx+1)];
        options = optimset('PlotFcns', {@optimplotx, @optimplotfval}, 'Display', 'iter', ...
        'FunValCheck', 'on');

        [epsilon_final, ~, ~, OUTPUT] = fzero(@(epsilon) optimal_rotation...
            (epsilon, C, L, Q, OG_nm, RD_nm, lambda, modal_amp_sim, t, init, 18000), epsilon_range, options);
        disp(OUTPUT);
    else
        disp('no sign flip detected');
        [~, idx] = min(abs(transfer));
        options = optimset('PlotFcns', {@optimplotx, @optimplotfval}, 'Display', 'iter', ...
        'FunValCheck', 'on');

        [epsilon_final, ~, ~, OUTPUT] = fminbnd(@(epsilon) abs(optimal_rotation...
            (epsilon, C, L, Q, OG_nm, RD_nm, lambda, modal_amp_sim, t, init, 18000)), epsilon(idx-1), epsilon(idx+1), options);
        disp(OUTPUT);
    end
    close all;

    % Final calculation of transformation matrix and new constant linear and
    % quadratic terms
    [~, X, C_til, L_til, Q_til] = ...
        optimal_rotation(epsilon_final, C, L, Q, OG_nm, RD_nm, lambda, modal_amp_sim, t, init, 18000);

    [pod_u_til, pod_v_til, pod_vor_til, modal_amp_raw_til] = ...
        basis_transform(pod_ut, pod_vt, pod_vort, modal_ampt, RD_nm, X);

    %% Temp
    Gal_coeff_til = [C_til L_til Q_til];
    Mod = true;

    Gal_coeff_til = ode_coefficients(RD_nm, Gal_coeff_til, Mod);
    options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

    tic1 = tic;
    [t, modal_amp_til] = ode113(@(t,y) system_odes_mod(t,y,Gal_coeff_til), t, ...
        modal_amp_sim(init,1:RD_nm), options);
    toc(tic1);
    %% Temp
    
    % TODO change
    if length(t) > 2
        sample_freq = 1/(t(2) - t(1));
        fprintf('Detected Sampling Frequency %6.2f Hz\n\n', sample_freq);
    else
        error('must provide tspan with a range');
    end

    % Prepare data
    plot_data.num_modes     = RD_nm;
    plot_data.direct        = direct;
    plot_data.init          = init;
    plot_data.pod_ut        = pod_u_til;
    plot_data.pod_vt        = pod_v_til;
    plot_data.pod_vort      = pod_vor_til;
    plot_data.dimensions    = dimensions;
    plot_data.fft_window    = fft_window;
    plot_data.u_scale       = u_scale;
    plot_data.l_scale       = l_scale;
    plot_data.plot_type     = plot_type;
    plot_data.sample_freq   = sample_freq;
    plot_data.id            = models;
    plot_data.modal_amp     = modal_amp_til;
    plot_data.t             = t;
    plot_data.x             = x;
    plot_data.y             = y;
    plot_data.bnd_idx       = bnd_idx;

    % Generate plots
    produce_plots(plot_data);

    results.X               = X;
    results.C_til           = C_til;
    results.L_til           = L_til;
    results.Q_til           = Q_til;
    results.pod_u_til       = pod_u_til;
    results.pod_v_til       = pod_v_til;
    results.modal_amp_til   = modal_amp_raw_til;
    results.epsilon_final   = epsilon_final;
    
    % Save important coefficients
    if save_mod == true
        save([direct '\Mod Galerkin Coeff\Coeff_' num2str(run_num) '_m' num2str(RD_nm) '_' models{i} '.mat'], ...
            'results');
    end
end

end



