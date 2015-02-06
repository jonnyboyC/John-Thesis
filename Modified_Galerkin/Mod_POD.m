function Mod_POD(varargin)
% TODO fill out help 

format long g
close all
clc;

%List of fields that will be checked
fields = {  'RD_nm',        'plot_type',    'save_mod', ...
            'init',         'line_range',   'direct' ,...
            'run_num',      'type',         'fft_window'};

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
type        = problem.type;
fft_window  = problem.fft_window;

% Check that parallel pool is ready
if isempty(gcp)
    parpool;
end

gcp();

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
pod_u   = vars.results.pod_u;           % streamwise pod modes
pod_v   = vars.results.pod_v;           % spanwise pod modes
u_scale     = vars.results.u_scale;     % velocity scaling
l_scale     = vars.results.l_scale;     % length scaling
lambda2     = vars.results.lambda2;     % eigenvalues of modes
modal_amp_raw = vars.results.modal_amp_raw;   % modal amplitude of each image in pod basis
dimensions  = vars.results.dimensions;  % dimensions of mesh
run_num     = vars.results.run_num;     % POD run numbers

clear vars

% Main Loop
for i = 1:size(type, 2)

    % TODO get load_data working 
    vars = load_data(type{i}, direct_Gal);
    
    % load variables for type
    
    q = vars.q;
    c = vars.c;
    l = vars.l;
    t = vars.t;
    OG_nm = vars.OG_nm;
    modal_amp = vars.modal_amp;

    clear vars

    % Truncate POD
    pod_ut = pod_u(:,1:OG_nm);
    pod_vt = pod_v(:,1:OG_nm);
    modal_amp_rawt = modal_amp_raw(:, 1:OG_nm);
    

    line_problem.c = c;
    line_problem.l = l;
    line_problem.q = q;
    line_problem.t = t;
    line_problem.init = init;
    line_problem.OG_nm = OG_nm;
    line_problem.RD_nm = RD_nm;
    line_problem.lambda2 = lambda2;
    line_problem.modal_amp = modal_amp;
    line_problem.line_range = line_range;

    [epsilon, transfer, flip_idx] = line_search(line_problem);

    if flip_idx ~= 0
        epsilon_range = [epsilon(flip_idx) epsilon(flip_idx+1)];
        options = optimset('PlotFcns', {@optimplotx, @optimplotfval}, 'Display', 'iter', ...
        'FunValCheck', 'on');

        [epsilon_final, ~, ~, OUTPUT] = fzero(@(epsilon) optimal_rotation...
            (epsilon, c, l, q, OG_nm, RD_nm, lambda2, modal_amp, t, init, 18000), epsilon_range, options);
        close all;
        disp(OUTPUT);
    else
        disp('no sign flip detected');
        [~, idx] = min(abs(transfer));
        options = optimset('PlotFcns', {@optimplotx, @optimplotfval}, 'Display', 'iter', ...
        'FunValCheck', 'on');

        [epsilon_final, ~, ~, OUTPUT] = fminbnd(@(epsilon) abs(optimal_rotation...
            (epsilon, c, l, q, OG_nm, RD_nm, lambda2, modal_amp, t, init, 18000)), epsilon(idx-1), epsilon(idx+1), options);
        disp(OUTPUT);
    end

    % Final calculation of transformation matrix and new constant linear and
    % quadratic terms
    [~, X, C_til, L_til, Q_til] = ...
        optimal_rotation(epsilon_final, c, l, q, OG_nm, RD_nm, lambda2, modal_amp, t, init, 18000);

    [pod_u_til, pod_v_til, modal_amp_raw_til] = ...
        basis_transform(pod_ut, pod_vt, modal_amp_rawt, RD_nm, X);

    %% Temp
    Gal_coeff_til = [C_til L_til Q_til];

    reduced_model_coeff_vis = -ode_coefficients(RD_nm, RD_nm, Gal_coeff_til);
    options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

    tic1 = tic;
    [t, modal_amp_til] = ode113(@(t,y) system_odes(t,y,reduced_model_coeff_vis), t, ...
        modal_amp(init,1:RD_nm), options);
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
    plot_data.dimensions    = dimensions;
    plot_data.fft_window    = fft_window;
    plot_data.u_scale       = u_scale;
    plot_data.l_scale       = l_scale;
    plot_data.plot_type     = plot_type;
    plot_data.sample_freq   = sample_freq;
    plot_data.id            = type;
    plot_data.modal_amp     = modal_amp_til;
    plot_data.t             = t;

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
        save([direct '\Mod Galerkin Coeff\Coeff_' num2str(run_num) '_m' num2str(RD_nm) '_' type{i} '.mat'], ...
            'results');
    end
end

end



