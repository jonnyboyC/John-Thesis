function problem = setdefaults_score(problem)

% Default for run_num
if isempty(problem.run_num) || (~isscalar(problem.run_num) && ~ischar(problem.run_num))
    fprintf('Using default value for run_num\nproblem.run_num = "first"\n\n');
    problem.run_num = 'first';        % use 10 modes
end

% Default for direct
if isempty(problem.direct) || ~ischar(problem.direct) 
    fprintf('Using default values for direct\nproblem.direct = ""\n\n');
    problem.direct = '';    % User prompt
end

% Default for num_pod 
if isempty(problem.num_clusters) || ~isnumeric(problem.num_clusters) && ~iscell(problem.num_clusters) 
    fprintf('Using default value for num_modesG\nproblem.num_modesG = 10\n\n');
    problem.num_clusters = 10;        % use 10 modes
end

% Default for num_cores
if isempty(problem.num_cores) || ~isscalar(problem.num_cores) || ...
        (ischar(problem.num_cores) && strcmp(problem.num_cores, 'auto'))
    fprintf('Using default value for num_cores\nproblem.num_cores = "num_cores"\n\n');
    problem.num_cores = 'auto';        % set to system max
end

end