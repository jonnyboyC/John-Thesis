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

% Prompt User for folder if directory is not provided
[direct_POD, direct] = prompt_folder('POD', run_num);

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

end