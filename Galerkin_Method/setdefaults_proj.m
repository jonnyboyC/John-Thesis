function problem = setdefaults_proj(problem)
% Function to set defaults for Galerkin_Proj

% Default for num_pod 
if isempty(problem.num_modesG) || ~isnumeric(problem.num_modesG) && ~iscell(problem.num_modesG) 
    fprintf('Using default value for num_modesG\nproblem.num_modesG = 10\n\n');
    problem.num_modesG = 10;        % use 10 modes
end

% Default for num_cores
if isempty(problem.num_cores) || ~isscalar(problem.num_cores) || ...
        (ischar(problem.num_cores) && strcmp(problem.num_cores, 'auto'))
    fprintf('Using default value for num_cores\nproblem.num_cores = "num_cores"\n\n');
    problem.num_cores = 'auto';        % set to system max
end

% Default for save_coef
if isempty(problem.score_model) || ~islogical(problem.score_model)
    fprintf('Using default value for score_model\nproblem.score_model = true\n\n');
    problem.score_model = true;        % save projection values
end

if iscell(problem.num_modesG)
    idx = cellfun(@isnumeric, problem.num_modesG);
    problem.num_modesG(~idx) = [];
end

% Default for int_time
if isempty(problem.int_time) || ~isscalar(problem.int_time)
    fprintf('Using default value for int_time\nproblem.int_time = 3600\n\n');
    problem.int_time = 3600;    
end

% Default for plot_type
if isempty(problem.plot_type) && ~iscell(problem.plot_type);
    fprintf('Using default value for plot_type\nproblem.plot_type = {"amp", "fft", "energy"}\n\n');
    problem.plot_type = {'amp', 'fft', 'energy'};      % plot types to be used
end

% Check to make sure incorrect strings are not passed
correct = {'fft', 'amp', 'energy', 'video', 'video stream'};
problem.plot_type = list_check(problem.plot_type, correct, 'plot_type');

% Default for save_coef
if isempty(problem.save_coef) || ~islogical(problem.save_coef)
    fprintf('Using default value for save_coef\nproblem.save_coef = true\n\n');
    problem.save_coef = true;        % save projection values
end

% Default for override_coef
if isempty(problem.override_coef) || ~islogical(problem.override_coef)
    fprintf('Using default value for override_coef\nproblem.override_coef = false\n\n');
    problem.override_coef = false;     % override previous coefficients
end

% Default for tspan
if isempty(problem.tspan) || ((~(isnumeric(problem.tspan) && ...
        size(problem.tspan, 1) == 1 && size(problem.tspan,2) > 1)) && ~iscell(problem.tspan))
    fprintf('Using default value for tspan\nproblem.tspan = 0:0.01:100\n\n');
    problem.tspan = {'test'};     % time range of integration
end

% Check to make sure incorrect strings are not passed
if iscell(problem.tspan)
    correct = {'test'};
    problem.tspan(1) = list_check(problem.tspan(1), correct, 'tspan');
end

if iscell(problem.tspan) && size(problem.tspan, 2) >= 2 && ~isnumeric(problem.tspan{2})
   fprintf('Using default value for tspan multiplier\nproblem.tspan{2} = 1\n\n');
   problem.tspan = problem.tspan(1);
end

if iscell(problem.tspan) && size(problem.tspan, 2) == 3 && ~isnumeric(problem.tspan{3})
   fprintf('Using default value for tspan extender\nproblem.tspan{2} = 1\n\n');
   problem.tspan = problem.tspan(1);
end

% Default for init
if isempty(problem.init) || ~isscalar(problem.init)
    fprintf('Using default value for init\nproblem.init = 1\n\n');
    problem.init = 1;   % initial condition snapshot
end

% Default for direct
if isempty(problem.direct) || ~ischar(problem.direct) 
    fprintf('Using default values for direct\nproblem.direct = ""\n\n');
    problem.direct = '';    % User prompt
end

% Default for Re0
if isempty(problem.Re0_gen) || (~isscalar(problem.Re0_gen) && ~isa(problem.Re0_gen, 'function_handle'))
    fprintf('Using default values for Re0_gen\nproblem.Re0_gen = @Re0_gen\n\n');
    problem.Re0_gen = @Re0_gen_shear;   	% Claimed high speed side
end

% Default for odesolver
if isempty(problem.odesolver) || (~isscalar(problem.odesolver) && ~isa(problem.odesolver, 'function_handle'))
    fprintf('Using default values for odesolver\nproblem.odesolver = @ode113\n\n');
    problem.odesolver = @ode113;   	% Selected ode solver method
end

% Default for fft_window
if isempty(problem.fft_window) || (isnumeric(problem.fft_window) && ...
        size(problem.fft_window, 1) == 1 && size(problem.fft_window,2) == 2)
    fprintf('Using default value for fft_window\nproblem.fft_window = [0 2000]\n\n');
    problem.fft_window = [0 2000];     % time range of integration
end

% Default for dissapation
if isempty(problem.dissapation) && ~iscell(problem.dissapation);
    fprintf('Using default value for dissapation\nproblem.dissapation = {"Least Squares", "Averaged"}\n\n');
    problem.dissapation = {'Least Squares', 'Averaged'};      % plot types to be used
end

% Check to make sure incorrect strings are not passed
correct = {'Least Squares', 'Averaged'};
problem.dissapation = list_check(problem.dissapation, correct, 'dissapation');

% Default for run_num
if isempty(problem.run_num) || (~isscalar(problem.run_num) && ~ischar(problem.run_num))
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

% Default for use time integration
if isempty(problem.time_int) || ~islogical(problem.time_int)
    fprintf('Using default value for time_int\nproblem.time_int = true\n\n');
    problem.time_int = true;      % time integrate galerkin systems
end

% Default for use coefficient calculation
if isempty(problem.calc_coef) || ~islogical(problem.calc_coef)
    fprintf('Using default value for calc_coef\nproblem.calc_coef = true\n\n');
    problem.calc_coef = true;      % time integrate galerkin systems
end

% Default for use chunks
if isempty(problem.use_chunks) || ~islogical(problem.use_chunks)
    fprintf('Using default value for use_chunks\nproblem.use_chunks = false\n\n');
    problem.use_chunks = false;     % use chunks to calculate quadractic terms
end

end

