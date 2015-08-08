function res_socres = score_all(varargin)
% SCORE_ALL perform generate and score surrogate markov models for each
% model present in a given folder

fields = {  'run_num',  'direct',   'num_clusters' , 'num_cores'};
    
% Parse problem structure provided to set it up correctly
if nargin == 1
    problem = parse_inputs(fields, @setdefaults_score, varargin{1});
else
    fprintf('Provide a single structure as input, use help Galerkin_Proj for information.\n');
    fprintf('Using Defaults\n\n');
    problem = parse_inputs(fields, @setdefaults_score);
end

% Create more readable names
run_num = problem.run_num;
direct  = problem.direct;
num_clusters = problem.num_clusters;
num_cores = problem.num_cores;

% Setup MATLAB to a max number of cores
setup_cores(num_cores);

% Prompt User for folder if directory is not provided
if strcmp(direct, '');
    [direct_POD, direct] = prompt_folder('POD', run_num);
else
    [direct_POD, direct] = prompt_folder('POD', run_num, direct);
end

% Check folders are up to most recent format
update_folders(direct);

% Load POD variables
vars = load(direct_POD);

run_num = vars.results_pod.run_num;
modal_amp = vars.results_pod.modal_amp;

if isfield(vars, 'results_clust');
    temp_clusters = vars.results_clust.num_clusters; % number of clustered used
    if temp_clusters == num_clusters
        km = vars.results_clust.km; % k-means cluster information
        gmm = vars.results_clust.gmm; % gaussian mixture model cluster information
    else
        % Fill with stub if dissimilar
        km = struct;
        gmm = struct;
    end
else
    % Fill with stub if empty
    km = struct;
    gmm = struct;
end

galerkin_folder = [filesep 'Galerkin Coeff' filesep];
galerkin_folder = [direct, galerkin_folder];
folders = dir(galerkin_folder);

for i = 3:length(folders)
    if ~folders.isdir
        continue;
    end
    
end


