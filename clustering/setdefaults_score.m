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

% Default for direct
if isempty(problem.score_mod) || ~ischar(problem.score_mod) 
    fprintf('Using default values for score_mod\nproblem.score_mod = ""\n\n');
    problem.score_mod = false;    % User prompt
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
end