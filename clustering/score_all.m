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

galerkin_path = [filesep 'Galerkin Coeff' filesep];
galerkin_path = [direct, galerkin_path];
files = dir(galerkin_path);

for i = 3:length(files)
    if ~files(i).isdir
        continue;
    end
    
    % If folder concat full path
    full_path = [galerkin_path, files(i).name];
    
    % Get first mat file in that folder
    file = get_file(run_num, full_path);
    
    vars = load([full_path filesep file]);
    
    
    % TODO have to get tspan from one dude that actually competed 
    % and backcalculate shit... god damn it
    
    % Ready score info Structure
    score_info.km = km;
    score_info.gmm = gmm;
    score_info.integration = vars.results_int.integration
    score_info.tspan = tspan;
    score_info.modal_amp = modal_amp;
    score_info.num_clusters = num_clusters;
    score_info.int_time = int_time;
    score_info.num_cores = num_cores;
    score_info.modes = modes;
    score_info.custom = custom;
    score_info.direct = direct;
    score_info.multiplier = multiplier;
    score_info.MOD = false;
    
    % score results
    [frob_km{i}, frob_gmm{i}, prob_km{i}, prob_gmm{i}, completed{i}] = ...
        score_model(score_info);
end


