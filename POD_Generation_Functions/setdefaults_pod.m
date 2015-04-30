function problem = setdefaults_pod(problem)
% Function to set defaults for POD_Gen

% Default for num_images 
if isempty(problem.num_images) || ~isscalar(problem.num_images)
    fprintf('Using default value for num_images\nproblem.num_images = 2000\n\n');
    problem.num_images = 2000;  % use 2000 images
end

% Default for load_raw
if isempty(problem.load_raw) || ~islogical(problem.load_raw)
    fprintf('Using default value for load_raw\nproblem.load_raw = true\n\n');
    problem.load_raw = true;        % read raw data
end

% Default for save_pod
if isempty(problem.save_pod) || ~islogical(problem.save_pod)
    fprintf('Using default value for save_pod\nproblem.save_pod = true\n\n');
    problem.save_pod = true;        % save pod results
end

% Default for save_figures
if  ~iscell(problem.save_figures)
    fprintf('Using default value for save_figures\nproblem.save_figures = {}\n\n');
    problem.save_figures = {};     % Save as figure
end

% Check to make sure incorrect strings are not passed
correct = {'fig', 'jpg'};
correct_members = ismember(problem.save_figures, correct);
for i = 1:size(correct_members,2)
    if ~correct_members(i)
        fprintf('%s is not a correct input\n', problem.save_figure{i});
    end
end
problem.save_figures = problem.save_figures(correct_members);

% Default for direct
if isempty(problem.direct) || ~ischar(problem.direct)
    fprintf('Using default values for direct\nproblem.direct = ""\n\n');
    problem.direct = '';
end

% Default for l_scale
if isempty(problem.l_scale) || ~isscalar(problem.l_scale)
    fprintf('Using default values for l_scale\nproblem.l_scale = 1\n\n');
    problem.l_scale = 1;    % Length Scale for problem
end

% Default for u_scale_gen
if isempty(problem.u_scale_gen) || (~isscalar(problem.u_scale_gen) && ~isa(problem.u_scale_gen, 'function_handle'))
    fprintf('Using default values for u_scale_gen\nproblem.u_scale_gen = 1\n\n');
    problem.u_scale_gen = 1;  	% Length scale for problem
end

% Default for flip
if isempty(problem.flip) || ~isequal(size(problem.flip), [1, 4]) || ~all(arrayfun(@islogical, problem.flip)) 
    fprintf('Using default value for flip\nproblem.flip = [false, false, false, false]\n\n');
    problem.flip = [false, false, false, false];        % save pod results
end

% Default for update_bnds
if isempty(problem.update_bnds) || ~islogical(problem.update_bnds)
    fprintf('Using default value for new_mask\nproblem.update_bnds = false\n\n');
    problem.update_bnds = false;        % save pod results
end

% Default for num_clusters
if isempty(problem.num_clusters) || ~isscalar(problem.num_clusters)
    fprintf('Using default values for num_clusters\nproblem.num_clusters = 10\n\n');
    problem.num_clusters = 10;  	% Length scale for problem
end

% Default for num_clusters
if isempty(problem.exp_sampling_rate) || ~isscalar(problem.exp_sampling_rate)
    fprintf('Using default values for exp_sampling_rate\nproblem.exp_sampling_rate = 5\n\n');
    problem.exp_sampling_rate = 5;  	% Length scale for problem
end


end

