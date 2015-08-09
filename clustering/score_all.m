function res_scores = score_all(varargin)
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
exp_sampling_rate = vars.results_pod.exp_sampling_rate;

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

results_scores = struct;
galerkin_path = [filesep 'Galerkin Coeff' filesep];
galerkin_path = [direct, galerkin_path];
files = dir(galerkin_path);

% Skip . and ..
% for i = 3:length(files)
%     close all
%     if ~files(i).isdir
%         continue;
%     end
%     
%     % If folder concat full path
%     full_path = [galerkin_path, files(i).name];
%     
%     % Get first mat file in that folder
%     file = get_file(run_num, full_path);
%     
%     if strcmp(file, 'not found')
%         fprintf('folder %s has a not matching run_num for %d\n', files(i).name, run_num);
%         continue;
%     end
%     
%     vars = load([full_path filesep file]);
%         
%     % TODO have to get tspan from one dude that actually competed 
%     % and backcalculate shit... god damn it
%     
%     integration = vars.results_int.integration;
%     
%     [tspan, multiplier] = back_calc_tspan(exp_sampling_rate, integration, modal_amp);
%     
%     % Ready score info Structure
%     score_info.km = km;
%     score_info.gmm = gmm;
%     score_info.integration = integration;
%     score_info.tspan = tspan;
%     score_info.modal_amp = modal_amp;
%     score_info.num_clusters = num_clusters;
%     score_info.int_time = 3600;
%     score_info.num_cores = num_cores;
%     score_info.modes = vars.results_coef.modes;
%     if strcmp(files(i).name, 'custom')
%         score_info.custom = true;
%     else
%         score_info.custom = false;
%     end
%     score_info.direct = direct;
%     score_info.multiplier = multiplier;
%     score_info.MOD = false;
%     
%     % score results
%     [frob_km, frob_gmm, prob_km, prob_gmm, completed] = ...
%         score_model(score_info);
%         
%     m = flow_comps(frob_km);
%     models = flow_ncomps(frob_km);
%     for j = 1:models
%         s = flow_comps(frob_km.(m{j}));
%         sub_models = flow_ncomps(frob_km.(m{j}));
%         for k = 1:sub_models
%             results_scores.(files(i).name).(m{j}).(s{k}).('frob_km') = frob_km.(m{j}).(s{k});
%             results_scores.(files(i).name).(m{j}).(s{k}).('frob_gmm') = frob_gmm.(m{j}).(s{k});
%             results_scores.(files(i).name).(m{j}).(s{k}).('prob_km') = prob_km.(m{j}).(s{k});
%             results_scores.(files(i).name).(m{j}).(s{k}).('prob_gmm') = prob_gmm.(m{j}).(s{k});
%             results_scores.(files(i).name).(m{j}).(s{k}).('completed') = completed.(m{j}).(s{k});
%         end
%     end
% end

mod_path = [filesep 'Mod Galerkin Coeff' filesep];
mod_path = [direct, mod_path];
files = dir(mod_path);

for i = 3:length(files)
    close all
    if ~files(i).isdir
        continue;
    end
    
    % If folder concat full path
    full_path = [mod_path, files(i).name];
    
    % Get first mat file in that folder
    file = get_file(run_num, full_path);
    
    if strcmp(file, 'not found')
        fprintf('folder %s has a not matching run_num for %d\n', files(i).name, run_num);
        continue;
    end
    
    vars = load([full_path filesep file]);
        
    % TODO have to get tspan from one dude that actually competed 
    % and backcalculate shit... god damn it
    
    integration = vars.results_mod_int.integration;
    
    [tspan, multiplier] = back_calc_tspan(exp_sampling_rate, integration, modal_amp);
    
    % Ready score info Structure
    score_info.km = km;
    score_info.gmm = gmm;
    score_info.integration = integration;
    score_info.tspan = tspan;
    score_info.modal_amp = modal_amp;
    score_info.num_clusters = num_clusters;
    score_info.int_time = 3600;
    score_info.num_cores = num_cores;
    score_info.modes = vars.results_coef.modes;
    if strcmp(files(i).name, 'custom')
        score_info.custom = true;
    else
        score_info.custom = false;
    end
    score_info.direct = direct;
    score_info.multiplier = multiplier;
    score_info.MOD = true;
    
    % score results
    [frob_km, frob_gmm, prob_km, prob_gmm, completed] = ...
        score_model(score_info);
        
    m = flow_comps(frob_km);
    models = flow_ncomps(frob_km);
    for j = 1:models
        s = flow_comps(frob_km.(m{j}));
        sub_models = flow_ncomps(frob_km.(m{j})); 
        s_new = cell(size(s));
        for k = 1:sub_models
            s_new{i} = ['mod_' s{i}];
        end
        for k = 1:sub_models
            results_scores.(files(i).name).(m{j}).(s_new{k}).('frob_km') = frob_km.(m{j}).(s{k});
            results_scores.(files(i).name).(m{j}).(s_new{k}).('frob_gmm') = frob_gmm.(m{j}).(s{k});
            results_scores.(files(i).name).(m{j}).(s_new{k}).('prob_km') = prob_km.(m{j}).(s{k});
            results_scores.(files(i).name).(m{j}).(s_new{k}).('prob_gmm') = prob_gmm.(m{j}).(s{k});
            results_scores.(files(i).name).(m{j}).(s_new{k}).('completed') = completed.(m{j}).(s{k});
        end
    end
end

if nargout == 1
    res_scores = results_scores;
end

