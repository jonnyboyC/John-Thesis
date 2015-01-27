function problem = setdefaults_proj(problem)
% Function to set defaults for Galerkin_Proj

% Default for num_pod 
if isempty(problem.num_pods) || ~isscalar(problem.num_pods)
    fprintf('Using default value for num_pods\nproblem.num_pods = 1000\nn');
    problem.num_pods = 10;        % use 10 modes
end

% Default for plot_pred
if isempty(problem.plot_pred) || ~ischar(problem.plot_pred) || iscell(problem.plot_pred);
    fprintf('Using default value for plot_pred\nproblem.plot_pred = {"amp", "fft"}\nn');
    problem.plot_pred = {'amp', 'fft'};      % plot types to be used
end

% Check to make sure incorrect strings are not passed
correct = {'fft', 'amp', 'video'};
correct_members = ismember(problem.plot_pred, correct);
for i = 1:size(correct_members,2)
    if ~correct_members(i)
        fprintf('%s is not a correct input\n', problem.save_figure{i});
    end
end
problem.plot_pred = problem.plot_pred(correct_members);

% Default for save_coef
if isempty(problem.save_coef) || ~islogical(problem.save_coef)
    fprintf('Using default value for save_coef\nproblem.save_coef = true\nn');
    problem.save_coef = true;        % save projection values
end

% Default for override_coef
if isempty(problem.override_coef) || ~iscell(problem.override_coef)
    fprintf('Using default value for override_coef\nproblem.override_coef = false\nn');
    problem.override_coef = false;     % override previous coefficients
end

% Default for tspan
if isempty(problem.tspan) || ~(isnumeric(problem.tspan) && ...
        size(problem.tspan, 1) == 1 && size(problem.tspan,2) > 1)
    fprintf('Using default value for tspan\nproblem.tspan = 0:0.01:100\nn');
    problem.tspan = 0:0.01:100;     % time range of integration
end

% Default for init
if isempty(problem.init) || ~isscalar(problem.init)
    fprintf('Using default value for init\nproblem.init = 1\nn');
    problem.init = 1;   % initial condition snapshot
end

% Default for direct
if isempty(problem.direct) || ~ischar(problem.direct) || isempty(problem.run_num)
    fprintf('Using default values for direct\nproblem.direct = ""\nn');
    problem.direct = '';    % User prompt
end

% TODO currently assuming this function does not require additional inputs
% which is most likely not true if we are generating Reynolds number
% Default for Re0
if isempty(problem.Re0_gen) || (~isscalar(problem.Re0_gen) && ~isa(problem.Re0_gen, 'function_handle'))
    fprintf('Using default values for Re0_gen\nproblem.Re0_gen = @Re0_gen\nn');
    problem.Re0_gen = @Re0_gen;   	% Claimed high speed side
end

% Default for fft_window
if isempty(problem.fft_window) || (isnumeric(problem.fft_window) && ...
        size(problem.fft_window, 1) == 1 && size(problem.fft_window,2) == 2)
    fprintf('Using default value for fft_window\nproblem.fft_window = [0 2000]\nn');
    problem.fft_window = [0 2000];     % time range of integration
end

% Default for run_num
if isempty(problem.run_num) || ~isscalar(problem.run_num)
    fprintf('Using default value for run_num\nproblem.num_pods = 0\nn');
    problem.num_pods = 10;        % use 10 modes
end

end

