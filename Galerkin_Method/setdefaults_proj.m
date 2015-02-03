function problem = setdefaults_proj(problem)
% Function to set defaults for Galerkin_Proj

% Default for num_pod 
if isempty(problem.num_pods) || ~isscalar(problem.num_pods)
    fprintf('Using default value for num_pods\nproblem.num_pods = 10\n\n');
    problem.num_pods = 10;        % use 10 modes
end

% Default for plot_type
if isempty(problem.plot_type) || ~ischar(problem.plot_type) || iscell(problem.plot_type);
    fprintf('Using default value for plot_type\nproblem.plot_type = {"amp", "fft"}\n\n');
    problem.plot_type = {'amp', 'fft'};      % plot types to be used
end

% Check to make sure incorrect strings are not passed
correct = {'fft', 'amp', 'video'};
correct_members = ismember(problem.plot_type, correct);
for i = 1:size(correct_members,2)
    if ~correct_members(i)
        fprintf('%s is not a correct input for problem.plot_type\n', problem.plot_type{i});
    end
end
problem.plot_type = problem.plot_type(correct_members);

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
if isempty(problem.tspan) || ~(isnumeric(problem.tspan) && ...
        size(problem.tspan, 1) == 1 && size(problem.tspan,2) > 1)
    fprintf('Using default value for tspan\nproblem.tspan = 0:0.01:100\n\n');
    problem.tspan = 0:0.01:100;     % time range of integration
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

% TODO currently assuming this function does not require additional inputs
% which is most likely not true if we are generating Reynolds number
% Default for Re0
if isempty(problem.Re0_gen) || (~isscalar(problem.Re0_gen) && ~isa(problem.Re0_gen, 'function_handle'))
    fprintf('Using default values for Re0_gen\nproblem.Re0_gen = @Re0_gen\n\n');
    problem.Re0_gen = @Re0_gen_shear;   	% Claimed high speed side
end

% Default for fft_window
if isempty(problem.fft_window) || (isnumeric(problem.fft_window) && ...
        size(problem.fft_window, 1) == 1 && size(problem.fft_window,2) == 2)
    fprintf('Using default value for fft_window\nproblem.fft_window = [0 2000]\n\n');
    problem.fft_window = [0 2000];     % time range of integration
end

% Default for run_num
if isempty(problem.run_num) || ~isscalar(problem.run_num) || ~ischar(problem.run_num)
    fprintf('Using default value for run_num\nproblem.run_num = "first"\n\n');
    problem.run_num = 'first';        % use 10 modes
end

% Check to make sure incorrect strings are not passed
if ischar(problem.run_num)
    correct = {'first'};
    correct_members = ismember(problem.type, correct);
    if ~correct_members 
        fprintf('%s is not a correct input for problem.run_num\n', problem.run_num);
        problem.run_num = 'first';
    end
end

end

