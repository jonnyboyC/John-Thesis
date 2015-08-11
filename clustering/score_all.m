function res_scores = score_all(varargin)
% SCORE_ALL perform generate and score surrogate markov models for each
% model present in a given folder

format long g
close all
clc;

fields = {  'run_num',      'direct',       'num_clusters' , ...
            'num_cores',    'outlier_mode', 'score_mod', ...
            'target_freq'};
    
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
score_mod = problem.score_mod;
target_freq = problem.target_freq;
outlier_mode = problem.outlier_mode;

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
fix_gmm(direct);

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

results_scores  = struct;
TKE             = struct;
frequency       = struct;
amplitude1      = struct;
amplitude2      = struct;

galerkin_path = [filesep 'Galerkin Coeff' filesep];
galerkin_path = [direct, galerkin_path];
files = dir(galerkin_path);

% Skip . and ..
for i = 3:length(files)
    close all
    if ~files(i).isdir
        continue;
    end
    
    % If folder concat full path
    full_path = [galerkin_path, files(i).name];
    
    % Get first mat file in that folder
    file = get_file(run_num, full_path);
    
    if strcmp(file, 'not found')
        fprintf('folder %s has a not matching run_num for %d\n', files(i).name, run_num);
        continue;
    end
    
    fprintf('Calculating scores for %s\n', files(i).name)
    
    % Load variables
    vars = load([full_path filesep file]);
       
    % Unpack variables
    integration = vars.results_int.integration;
    modes = vars.results_coef.modes;
    MOD = false;
    sample_freq = vars.results_coef.sample_freq;
    
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
    score_info.outlier_mode = outlier_mode;
    score_info.modes = modes;
    score_info.direct = direct;
    score_info.multiplier = multiplier;
    score_info.MOD = MOD;
    
    % Can eventually remove this as well
    if strcmp(files(i).name, 'custom')
        score_info.custom = true;
    else
        score_info.custom = false;
    end
    
    % score results
    [frob_km, frob_gmm, like_km, like_gmm, completed] = ...
        score_model(score_info);
    
    % Get model names
    m = flow_comps(frob_km);
    models = flow_ncomps(frob_km);
    
    for j = 1:models
        % Get submodel names
        s = flow_comps(frob_km.(m{j}));
        sub_models = flow_ncomps(frob_km.(m{j}));
        
        for k = 1:sub_models  
            % Calculate TKE and pack results
            modal_amp_sim = integration.modal_amp.(m{j}).(s{k})(:,2:end);
            if target_freq ~= 0
                frequency = compare_freq(frequency, modal_amp_sim, sample_freq, ...
                    completed, target_freq, m{j}, s{k}, files(i));
            end
            
            % Caculate TKE and pack results
            TKE = calc_energy(TKE, modal_amp_sim, completed, modal_amp, m{j}, ...
                s{k}, files(i), modes);
            
            % Caculate values for first 2 modes
            [amplitude1, amplitude2] = compare_amp(amplitude1, amplitude2, ...
                modal_amp_sim, modal_amp, m{j}, s{k}, files(i), modes, completed);
            
            % Pack results of scores
            results_scores = pack_results(results_scores, completed, frob_km, ...
                frob_gmm, like_km, like_gmm, m{j}, s{k}, files(i));
        end
    end
end

mod_path = [filesep 'Mod Galerkin Coeff' filesep];
mod_path = [direct, mod_path];
files = dir(mod_path);

if score_mod
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
    
    % Load variables
    vars = load([full_path filesep file]);
    
    % Unpack variables
    integration = vars.results_mod_int.integration;
    modes = vars.results_mod_coef.modes;
    MOD = true;
    
    % Back calculate tspan and multiplier eventually can remove
    [tspan, multiplier] = back_calc_tspan(exp_sampling_rate, integration, modal_amp);
    
    % Ready score info structure
    score_info.km = km;
    score_info.gmm = gmm;
    score_info.integration = integration;
    score_info.tspan = tspan;
    score_info.modal_amp = modal_amp;
    score_info.num_clusters = num_clusters;
    score_info.int_time = 3600;
    score_info.num_cores = num_cores;
    score_info.outlier_mode = outlier_mode;
    score_info.modes = modes;
    score_info.modes = vars.results_mod_coef.modes;
    score_info.direct = direct;
    score_info.multiplier = multiplier;
    score_info.MOD = MOD;
    score_info.modal_amp_til = vars.results_mod_coef.system.modal_amp_til;
    
    % Can eventually remove this as well
    if strcmp(files(i).name, 'custom')
        score_info.custom = true;
    else
        score_info.custom = false;
    end
    
    % Score results
    [frob_km, frob_gmm, like_km, like_gmm, completed] = ...
        score_model(score_info);
        
    % Get Model names
    m = flow_comps(frob_km);
    models = flow_ncomps(frob_km);
    
    for j = 1:models
        % Get submodel names
        s = flow_comps(frob_km.(m{j}));
        sub_models = flow_ncomps(frob_km.(m{j})); 
        
        for k = 1:sub_models
            modal_amp_sim = integration.modal_amp.(m{j}).(s{k});
            if target_freq ~= 0
                frequency = compare_freq(frequency, modal_amp_sim, sample_freq, ...
                    completed, target_freq, m{j}, s{k}, files(i));
            end
            
            % Caculate TKE and pack results
            TKE = calc_energy(TKE, modal_amp_sim, completed, modal_amp, m{j}, ...
                s{k}, files(i), modes);
            
            % Caculate values for first 2 modes
            [amplitude1, amplitude2] = compare_amp(amplitude1, amplitude2, ...
                modal_amp_sim, modal_amp, m{j}, s{k}, files(i), modes, completed);
            
            % Pack results of scores
            results_scores = pack_results(results_scores, completed, frob_km, ...
                frob_gmm, like_km, like_gmm, m{j}, s{k}, files(i));
        end
    end
end
end

[fkm_list, model]  = merge_struct(results_scores, 'frob_km');
[fgmm_list, ~] = merge_struct(results_scores, 'frob_gmm');
[pkm_list, ~] = merge_struct(results_scores, 'like_km');
[pgmm_list, ~] = merge_struct(results_scores, 'like_gmm');
[tke_list, ~] = merge_struct(TKE, 'mean_diff');
[tke2_list, ~] = merge_struct(TKE, 'std_diff');
[fft_list, ~] = merge_struct(frequency, 'diff');
[fft2_list, ~] = merge_struct(frequency, 'mode');

% derp = strncmp(model, 'GM3', 3)

if nargout == 1
    res_scores = results_scores;
end

% TKE = sum(1/2*modal_amp'.^2);

end

