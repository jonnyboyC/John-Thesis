function problem = setdefaults_mod(problem)

% Default for Reduced num_modes 
if isempty(problem.RD_nm) || ~isscalar(problem.RD_nm)
    fprintf('Using default value for RD_nm\nproblem.RD_nm = 8\n\n');
    problem.RD_nm = 8;        % use 10 modes
end

% Default for type
if isempty(problem.plot_type) || ~ischar(problem.plot_type) || iscell(problem.plot_type);
    fprintf('Using default value for plot_type\nproblem.plot_type = {"amp", "fft"}\n\n');
    problem.plot_type = {'amp', 'fft'};      % plot_types to be used
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
    problem.line_range = 100;  % use 10 modes
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
if isempty(problem.type) || ~ischar(problem.type)
    fprintf('Using default value for type\nproblem.type = "vis1"\n\n');
    problem.type = 'vis1';      % previous galerkin type
end

% Check to make sure incorrect strings are not passed
correct = {'og', 'vis1', 'vis2'};
correct_members = ismember(problem.type, correct);
if ~correct_members
    fprintf('%s is not a correct input for problem.type\n', problem.type{i});
    problem.type = 'vis1';
end

end