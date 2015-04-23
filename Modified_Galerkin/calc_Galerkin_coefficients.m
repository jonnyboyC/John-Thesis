function calc_Galerkin_coefficients(varargin)

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
num_modesG      = problem.num_modesG+1; % Add mean flow
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
lambda      = vars.results.lambda;      % eigenvalues of modes
dimensions  = vars.results.dimensions;  % dimensions of mesh
vol_frac    = vars.results.vol_frac;    % mesh area size
bnd_idx     = vars.results.bnd_idx;     % location of boundaries
bnd_x       = vars.results.bnd_x;       % location of flow boundaries normal to x
bnd_y       = vars.results.bnd_y;       % location of flow boundaries normal to y
run_num     = vars.results.run_num;     % POD run numbers
cutoff      = vars.results.cutoff;      % number of modes at cutoff
modal_amp   = vars.results.modal_amp;   % modal amplitude  from raw data

% Get Reynolds number 
Re0 = Re0_gen(direct, u_scale, l_scale);  

% Kinematic visocity Assume values have been non-dimensionalized
vis = 1/Re0;

% Determine sampling frequency from provided tspan
if length(tspan) > 2
    sample_freq = 1/(tspan(2) - tspan(1));
    fprintf('Detected Sampling Frequency %6.2f Hz\n\n', sample_freq);
else
    error('must provide tspan with a range');
end

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
coef_problem.run_num        = run_num;
coef_problem.override_coef  = override_coef;
coef_problem.direct         = direct;

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

for i = 1:length(num_modesG)
    % Pull current number of modes to be investigated
    num_modes = num_modesG(i);

    % Created truncated pod basis
    pod_ut  = pod_u(:,1:num_modes);
    pod_vt  = pod_v(:,1:num_modes);

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
    
        % calculate coefficients detailed by Couplet
    if ismember('Couplet', dissapation);
        eddy{3,1,i} = viscous_dis_couplet(modal_amp, num_modes, lc{1,1}, qc{1,1}, vis);
        eddy{3,2,i} = 'Modal Base Couplet';

        eddy{4,1,i} = viscous_dis_couplet(modal_amp, num_modes, lc{2,1}, qc{2,1}, vis);
        eddy{4,2,i} = 'Modal Weak Couplet';
    end

    disp(eddy{3,1,i});
    
    % Calculate coefficeints detailed by Noack
    [eddy{5,1,i}, eddy{6,1,i}] = viscous_dis(modal_amp, num_modes, lambda, l{1,1,i}, q{1,1,i}, vis);
    eddy{5,2,i} = 'Modal Base Noack';
    eddy{6,2,i} = 'Global Base Noack';

    [eddy{7,1,i}, eddy{8,1,i}] = viscous_dis(modal_amp, num_modes, lambda, l{2,1,i}, q{2,1,i}, vis);
    eddy{7,2,i} = 'Modal Weak Noack';
    eddy{8,2,i} = 'Global Weak Noack';
    
    % duplicate TODO name this better
    eddy(linear_models+1:total_models,1,i) = eddy(base_models+1:linear_models,1,i);
    for j = base_models+1:linear_models
        eddy{j+linear_models-base_models,2,i} = ['NL ' eddy{j,2,i}];
    end
    
    % Prepare data
    results_coef.run_num = run_num;
    results_coef.l = squeeze(l(:,:,i));
    results_coef.q = squeeze(q(:,:,i));
    results_coef.eddy = squeeze(eddy(:,:,i));
    results_coef.vis = vis;
    results_coef.num_modesG = num_modes-1;
      results_coef.sample_freq = sample_freq;
    results_coef.linear_models = linear_models;
    results_coef.total_models = total_models;
    
    % Save relavent coefficients
    if save_coef == true
        if i ~= 1
            wait(futures)
        end
        pool = gcp;
        futures = parfeval(pool, @save_results, 0, results_coef, direct);
    end
end
end

function save_results(results_coef, direct)
    if ~exist([direct filesep 'Galerkin Coeff' filesep 'modes_' num2str(results.num_modesG)], 'dir') 
        mkdir([direct filesep 'Galerkin Coeff' filesep 'modes_' num2str(results.num_modesG)]);
    end
    save([direct filesep 'Galerkin Coeff' filesep 'modes_' num2str(results_coef.num_modesG-1) filesep 'Coefficients_run_'...
    num2str(results.run_num) '.mat'], 'results_coef', '-v7.3');
end