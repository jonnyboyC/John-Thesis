function problem = setdefaults_mod(problem)
% SETDEFAULTS_MOD setup the defaults for MOD_POD
%
%   problem = SETDEFAULTS_MOD(problem)

% Default for Reduced num_modes 
if isempty(problem.num_modes) || ~isscalar(problem.num_modes)
    fprintf('Using default value for num_modes\nproblem.num_modes = 8\n\n');
    problem.num_modes = 10;        % use 10 modes
end

% Default for type
if isempty(problem.plot_type) || ~ischar(problem.plot_type) || iscell(problem.plot_type);
    fprintf('Using default value for plot_type\nproblem.plot_type = {"amp", "fft", "energy"}\n\n');
    problem.plot_type = {'amp', 'fft', 'energy'};      % plot_types to be used
end

% Default for score
if isempty(problem.score) || ~islogical(problem.score)
    fprintf('Using default value for save_coef\nproblem.score = true\n\n');
    problem.score = true;        % save projection values
end

% Default for num_cores
if isempty(problem.num_cores) || ~isscalar(problem.num_cores) || ...
        (ischar(problem.num_cores) && strcmp(problem.num_cores, 'auto'))
    fprintf('Using default value for num_cores\nproblem.num_cores = "num_cores"\n\n');
    problem.num_cores = 'auto';        % set to system max
end

% Default for outlier_mode
if isempty(problem.outlier_mode) || ~isstruct(problem.outlier_mode)
    fprintf(['Using default value for outlier_mode\nproblem.outlier_mode.km = []' ...
       '\nproblem.outlier_mode.gm = [] \n\n']);
    problem.outlier_mode.km = [];
    problem.outlier_mode.gmm = [];
end

if isstruct(problem.outlier_mode)
    match = {'km', 'gmm'};
    fields =  fieldnames(problem.outlier_mode);
    rm_fields = fields(~ismember(fields, match));
    if ~isempty(rm_fields)
        for i = 1:size(rm_fields,1);
            fprintf('Removing field "%s"\n', rm_fields{i});
        end
        problem.outlier_mode = rmfield(problem.outlier_mode, rm_fields);
    end
    
    fields =  fieldnames(problem.outlier_mode);
    for i = 1:length(fields)
        if ~isscalar(problem.outlier_mode.(fields{i}))
           problem.outlier_mode.fields{i} = []; 
        end
    end
end

% Check to make sure incorrect strings are not passed
correct = {'fft', 'amp', 'energy', 'video', 'video stream'};
problem.plot_type = list_check(problem.plot_type, correct, 'plot_type');

% Default for save_coef
if isempty(problem.custom) || ~islogical(problem.custom)
    fprintf('Using default value for custom\nproblem.custom = false\n\n');
    problem.custom = false;        % save projection values
end

% Default for fft_window
if isempty(problem.fft_window) || (isnumeric(problem.fft_window) && ...
        size(problem.fft_window, 1) == 1 && size(problem.fft_window,2) == 2)
    fprintf('Using default value for fft_window\nproblem.fft_window = [0 2000]\n\n');
    problem.fft_window = [0 2000];     % time range of integration
end

% Default for tspan
if isempty(problem.tspan) || ((~(isnumeric(problem.tspan) && ...
        size(problem.tspan, 1) == 1 && size(problem.tspan,2) > 1)) && ~iscell(problem.tspan))
    fprintf('Using default value for tspan\nproblem.tspan = 0:0.01:100\n\n');
    problem.tspan = 0:0.0001:1;     % time range of integration
end

% Check to make sure incorrect strings are not passed
if iscell(problem.tspan)
    correct = {'test'};
    problem.tspan(1) = list_check({problem.tspan{1}}, correct, 'tspan');
end

if iscell(problem.tspan) && size(problem.tspan, 2) == 2 && ~isnumeric(problem.tspan{2})
   fprintf('Using default value for tspan multiplier\nproblem.tspan{2} = 1\n\n');
   problem.tspan = problem.tspan{1};
end

% Default for save_mod
if isempty(problem.save_mod) || ~islogical(problem.save_mod)
    fprintf('Using default value for save_mod\nproblem.save_mod = true\n\n');
    problem.save_mod = true;        % save modified projection values
end

% Default for init
if isempty(problem.init) || ~isscalar(problem.init)
    fprintf('Using default value for init\nproblem.init = 1\n\n');
    problem.init = 1;   % initial condition snapshot
end

% Default for line range 
if isempty(problem.line_range) || ~isscalar(problem.line_range)
    fprintf('Using default value for line_range\nproblem.line_range = 100\n\n');
    problem.line_range = 200;  % use 10 modes
end

% Default for direct
if isempty(problem.direct) || ~ischar(problem.direct) 
    fprintf('Using default values for direct\nproblem.direct = ""\n\n');
    problem.direct = '';    % User prompt
end

% Default for run_num
if isempty(problem.run_num) || ~isscalar(problem.run_num) || ~ischar(problem.run_num)
    fprintf('Using default value for run_num\nproblem.run_num = "first"\n\n');
    problem.run_num = 'first';        % use 10 modes
end

% Check to make sure incorrect strings are not passed
if ischar(problem.run_num)
    correct = {'first'};
    correct_members = ismember(problem.run_num, correct);
    if ~correct_members 
        fprintf('%s is not a correct input for problem.run_num\n', problem.run_num);
        problem.run_num = 'first';
    end
end

% Default for previous galerkin type
if isempty(problem.models) || (~iscell(problem.models) && ~ischar(problem.models))
    fprintf('Using default value for models\nproblem.models = {"GM"}\n\n');
    problem.models = {'GM'};            % Generating basis model
end

if iscell(problem.models)
    % Check to make sure incorrect strings are not passed
    correct = {'GM', 'GM1', 'GM2', 'GM3', 'all'};
    problem.models = list_check(problem.models, correct, 'models');
end

% Default for previous galerkin type
if isempty(problem.submodels) || ~iscell(problem.submodels) || ~ischar(problem.submodels)
    fprintf('Using default value for submodels\nproblem.submodels = {"base", "weak"}\n\n');
    problem.submodels = {'base', 'weak'}; % Generating basis submodel
end

if iscell(problem.submodels)
    % Check to make sure incorrect strings are not passed
    correct = {'base', 'weak', 'all'};
    problem.submodels = list_check(problem.submodels, correct, 'models');
end

% Check to make sure incorrect strings are not passed
if ischar(problem.submodels)
    correct = {'all'};
    correct_members = ismember(problem.submodels, correct);
    if ~correct_members 
        fprintf('%s is not a correct input for problem.submodels\n', problem.submodels);
        problem.submodels = 'all';
    end
end

% Default of basis_modes
if isempty(problem.basis_modes) || ~isscalar(problem.basis_modes) || ~ischar(problem.basis_modes)
    fprintf('Using default value for basis_modes\nproblem.basis_modes = "double"\n\n');
    problem.basis_modes = 2*problem.num_modes;
end