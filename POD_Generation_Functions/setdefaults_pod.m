function problem = setdefaults_pod(problem)
% Function to set defaults for POD_Gen

% Default for num_images 
if isempty(problem.num_images) || ~isscalar(problem.num_images)
    fprintf('Using default value for num_images\nproblem.num_images = 2000\n\n');
    problem.num_images = 2000;          % use 2000 images
end

% Default for load_raw
if isempty(problem.load_raw) || ~islogical(problem.load_raw)
    fprintf('Using default value for load_raw\nproblem.load_raw = true\n\n');
    problem.load_raw = true;            % read raw data
end

% Default for load_raw
if isempty(problem.streamlines) || ~islogical(problem.streamlines)
    fprintf('Using default value for streamlines\nproblem.streamlines = false\n\n');
    problem.streamlines = true;        % plot pod modes with quivers
end

% Default for num_cores
if isempty(problem.num_cores) || ~isscalar(problem.num_cores) || ...
        (ischar(problem.num_cores) && strcmp(problem.num_cores, 'auto'))
    fprintf('Using default value for num_cores\nproblem.num_cores = "num_cores"\n\n');
    problem.num_cores = 'auto';        % set to system max
end

% Default for grid_direct
if isempty(problem.grid_direct) || ~iscell(problem.grid_direct)
    fprintf('Using default value for grid_direct\nproblem.grid_direct = {}\n\n');
    problem.grid_direct = {};       % TODO
end

% Default for save_pod
if isempty(problem.save_pod) || ~islogical(problem.save_pod)
    fprintf('Using default value for save_pod\nproblem.save_pod = true\n\n');
    problem.save_pod = true;            % save pod results
end

% Default for cluster
if isempty(problem.cluster) || ~islogical(problem.cluster)
    fprintf('Using default value for cluster\nproblem.cluster = true\n\n');
    problem.cluster = true;             % cluster modal amplitudes
end

% Default for filter
if isempty(problem.filter) || ~islogical(problem.filter)
    fprintf('Using default value for cluster\nproblem.filter = false\n\n');
    problem.filter = false;             % filter using a guide filter
end

% Default for average_mesh
if isempty(problem.average_mesh) || ~islogical(problem.average_mesh)
    fprintf('Using default value for cluster\nproblem.average_mesh = false\n\n');
    problem.average_mesh = false;       % condense mesh using 2x2 pixels
end

% Default for non_dim
if isempty(problem.non_dim) || ~islogical(problem.non_dim)
    fprintf('Using default value for cluster\nproblem.non_dim = false\n\n');
    problem.non_dim = false;            % non-dimensionalize
end

% Default for image_range
if ~iscell(problem.image_range) || any(~cellfun(@isvector, problem.image_range))
     fprintf('Using default value for image_range\nproblem.image_range = {}\n\n');
     problem.image_range = {};     % perform no cropping of image
end

% Default for save_figures
if ~iscell(problem.save_figures)
    fprintf('Using default value for save_figures\nproblem.save_figures = {"fig", "png"}\n\n');
    problem.save_figures = {'fig', 'png'};     % Save figures as .fig and .png
end

% Check to make sure incorrect strings are not passed
correct = {'fig', 'jpg', 'png'};
problem.save_figures = list_check(problem.save_figures, correct, 'save_figures');

% Default for xy_units
if isempty(problem.xy_units) || ~ischar(problem.xy_units)
    fprintf('Using default values for direct\nproblem.xy_units = "mm"\n\n');
    problem.xy_units = 'mm';            % assume data is in millimeters
end

% Correct if invalid unit is given
correct = {'mm', 'm'};
found = strcmp(problem.xy_units, correct);
if isempty(found)
    problem.xy_units = 'mm';
end

% Default for direct
if isempty(problem.direct) || ~ischar(problem.direct)
    fprintf('Using default values for direct\nproblem.direct = ""\n\n');
    problem.direct = '';                % prompt user for directory
end

% Default for l_scale
if isempty(problem.l_scale) || ~isscalar(problem.l_scale)
    fprintf('Using default values for l_scale\nproblem.l_scale = 1\n\n');
    problem.l_scale = 1;                % assume characteristic length of 1
end

% Default for load_handle
if isempty(problem.load_handle) || ~isa(problem.load_handle, 'function_handle')
    fprintf('Using default values for load_handle\nproblem.load_handle = @load_LaVision\n\n');
    problem.load_handle = @load_LaVision;   % Assume the more recent laVision format
end

% Default for u_scale_gen
if isempty(problem.u_scale_gen) || (~isscalar(problem.u_scale_gen) && ~isa(problem.u_scale_gen, 'function_handle'))
    fprintf('Using default values for u_scale_gen\nproblem.u_scale_gen = 1\n\n');
    problem.u_scale_gen = 1;            % assume characteristic velocity of 1
end

% Default for flow_flip
if  ~iscell(problem.flow_flip)
    fprintf('Using default value for flow_flip\nproblem.flow_flip = {}\n\n');
    problem.flow_flip = {};     % Save figures as .fig and .png
end

% Check entries are correct
correct = {'x', 'y', 'z', 'u', 'v', 'w'};
problem.flow_flip = list_check(problem.flow_flip, correct, 'flow_flip');

% Default for update_bnds
if isempty(problem.update_bnds) || ~islogical(problem.update_bnds)
    fprintf('Using default value for new_mask\nproblem.update_bnds = false\n\n');
    problem.update_bnds = false;        % update flow boundaries
end

% Default for num_clusters
if isempty(problem.num_clusters) || ~isscalar(problem.num_clusters)
    fprintf('Using default values for num_clusters\nproblem.num_clusters = 10\n\n');
    problem.num_clusters = 10;  	% Length scale for problem
end

% Default for experimental sampling rate
if isempty(problem.exp_sampling_rate) || ~isscalar(problem.exp_sampling_rate)
    fprintf('Using default values for exp_sampling_rate\nproblem.exp_sampling_rate = 5\n\n');
    problem.exp_sampling_rate = 5;  	% Length scale for problem
end


% Default for load_only
if isempty(problem.load_only) || ~islogical(problem.load_only)
    fprintf('Using default value for read_only\nproblem.load_only = false\n\n');
    problem.load_only = false;      % load to produce mask and mat file
end

% Override previous flags load only is set
if problem.load_only == true
    if problem.load_raw == false
        problem.load_raw = true;
        fprintf(['Overriding setting for load_raw because load_only = true' ...
            '\nproblem.load_raw = true\n\n']); 
    end
    if problem.update_bnds == false
        problem.update_bnds = true;
        fprintf(['Overriding setting for update_bnds because load_only = true' ...
            '\nproblem.update_bnds = true\n\n']);
    end
end
end

