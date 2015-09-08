function res_scores = Score_all(varargin)
% SCORE_ALL  Produces scores for all modes of a given data set, currently
% genrates cluster at runtime so scores may vary slightly between runs
%
%   SCORE_ALL() Generates scores for all models of a given data set prompts
%   user for analysis folder and using default values detailed below
%
%   SCORE_ALL(problem) Gerneate scores for all models with flags specified 
%   by the structure problem specified below
%
%   res_scores = SCORE_ALL(problem) return a structure containing
%   information about for each model
%
%   problem.run_num = 'first
%   Specify which run this Galerkn Projection should be based from, default
%   is to use the most recent
%  
%   problem.direct = ''
%   Specify directory that will be searched for POD data, default is to
%   prompt user
%
%   problem.num_clusters = 10
%   Set the number of clsuter that should be determined.
% 
%   problem.num_cores  = 'auto'
%   Set the max number of cores to be used for POD_Gen, auto will set this
%   to the number of cores in the computer.
%
%   problem.outlier_mode.km = []
%   problem.outlier_mode.gmm = []
%   Set the distance of the outlier mode for km and gmm clustering, default
%   is to include no outlier mode, see John Chabot Thesis 2015 for details
%
%   problem.score_mod = false
%   Boolean value to determine if basis transformed models will be scored
%
%   problem.target_freq = 0 
%   Set a target frequency to see how close a model gets to a particular
%   frequency value, default is to not check
%
%   problem.phase = []
%   Currently somewhat of a stub, default is to not check phase
%
%   problem.steady = false
%   Boolean flag to restrict models of interest to just those that do not
%   come to a fixed point, default is to include all models.

format long g
close all
clc;

% Problem structure fields
fields = {  'run_num',      'direct',       'num_clusters' , ...
            'num_cores',    'outlier_mode', 'score_mod', ...
            'target_freq',  'phase',        'steady'};
    
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
phase = problem.phase;
steady = problem.steady;

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
    
    %%%%%%%%%%%%%%%%%TODO reverse post thesis%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     temp_clusters = vars.results_clust.num_clusters; % number of clustered used
%     if temp_clusters == num_clusters
%         km = vars.results_clust.km; % k-means cluster information
%         gmm = vars.results_clust.gmm; % gaussian mixture model cluster information
%     else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Fill with stub if dissimilar
        km = struct;
        gmm = struct;
%     end
else
    % Fill with stub if empty
    km = struct;
    gmm = struct;
end

results_scores  = struct;
TKE             = struct;
amplitude1      = struct;
amplitude2      = struct;
phase_comp      = struct;
frequency       = struct;


galerkin_path = [filesep 'Galerkin Coeff' filesep];
galerkin_path = [direct, galerkin_path];
files = dir(galerkin_path);

% % Skip . and ..
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
    sample_freq = 8000;%vars.results_coef.sample_freq;
    
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
    [frob_km, frob_gmm, like_km, like_gmm, completed, km_steady, gmm_steady] = ...
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
            
            % Continue if we find an error
            if isstruct(modal_amp_sim)
                continue
            end
            
            % If requested test frequency
            if target_freq ~= 0
                frequency = compare_freq(frequency, modal_amp_sim, sample_freq, ...
                    target_freq, completed, km_steady, gmm_steady, m{j}, s{k}, files(i));
            end
            
            % If request test phase
            if ~isempty(phase)
                phase_comp = phase_comparision(phase_comp, phase, modal_amp_sim, ...
                    modal_amp, modes, completed, km_steady, gmm_steady, m{j}, s{k}, files(i));
            end
            
            % Caculate TKE and pack results
            TKE = calc_energy(TKE, modal_amp_sim, modal_amp,modes ,completed, ...
                km_steady, gmm_steady, m{j}, s{k}, files(i));
            
            % Caculate values for first 2 modes
            [amplitude1, amplitude2] = compare_amp(amplitude1, amplitude2, ...
                modal_amp_sim, modal_amp, modes, completed, km_steady, gmm_steady, m{j}, s{k}, files(i));
            
            % Pack results of scores
            results_scores = pack_results(results_scores, frob_km, frob_gmm, ...
                like_km, like_gmm,  completed, km_steady, gmm_steady, m{j}, s{k}, files(i));
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
    if isfield(vars.results_mod_int, 'integration') || ...
            isfield(vars.results_mod_coef, 'modes') || ...
            isfield(vars.results_mod_coef, 'system');
        integration = vars.results_mod_int.integration;
        modes = vars.results_mod_coef.modes;
        modal_amp_til = vars.results_mod_coef.system.modal_amp_til;
    else
        continue;
    end
    
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
    score_info.direct = direct;
    score_info.multiplier = multiplier;
    score_info.MOD = MOD;
    score_info.modal_amp_til = modal_amp_til;
    
    % Can eventually remove this as well
    if strcmp(files(i).name, 'custom')
        score_info.custom = true;
    else
        score_info.custom = false;
    end
    
    % score results
    [frob_km, frob_gmm, like_km, like_gmm, completed, km_steady, gmm_steady] = ...
        score_model(score_info);
        
    % Get Model names
    m = flow_comps(frob_km);
    models = flow_ncomps(frob_km);
    
    for j = 1:models
        % Get submodel names
        s = flow_comps(frob_km.(m{j}));
        sub_models = flow_ncomps(frob_km.(m{j})); 
        
        for k = 1:sub_models
            % Calculate TKE and pack results
            modal_amp_sim = integration.modal_amp.(m{j}).(s{k})(:,2:end);
            
            % Continue if we find an error
            if isstruct(modal_amp_sim)
                continue
            end
            
            % If requested test frequency
            if target_freq ~= 0
                frequency = compare_freq(frequency, modal_amp_sim, sample_freq, ...
                    target_freq, completed, km_steady, gmm_steady, m{j}, s{k}, files(i));
            end
            
            % Caculate TKE and pack results
            TKE = calc_energy(TKE, modal_amp_sim, modal_amp,modes ,completed, ...
                km_steady, gmm_steady, m{j}, s{k}, files(i));
            
            % Caculate values for first 2 modes
            [amplitude1, amplitude2] = compare_amp(amplitude1, amplitude2, ...
                modal_amp_sim, modal_amp, modes, completed, km_steady, gmm_steady, m{j}, s{k}, files(i));
            
            % Pack results of scores
            results_scores = pack_results(results_scores, frob_km, frob_gmm, ...
                like_km, like_gmm,  completed, km_steady, gmm_steady, m{j}, s{k}, files(i));
        end
    end
end
end


clusters_info.results_scores    = results_scores;
clusters_info.TKE               = TKE;
clusters_info.frequency         = frequency;
clusters_info.phase_comp        = phase_comp;
clusters_info.amplitude1        = amplitude1;
clusters_info.amplitude2        = amplitude2;
clusters_info.num_clusters      = num_clusters;
clusters_info.direct            = direct;
clusters_info.steady            = steady;

relations = cluster_plots(clusters_info);

results_scores.relations = relations;
if steady
    direct_ext = [direct filesep 'Scores' filesep 'steady' filesep num2str(num_clusters) 'cluster'];
else
    direct_ext = [direct filesep 'Scores' filesep 'not steady' filesep num2str(num_clusters) 'cluster'];
end

if ~exist(direct_ext, 'dir') 
    mkdir(direct_ext);
end

save([direct_ext filesep 'scores_' num2str(run_num)], 'results_scores', '-v7.3');

if nargout == 1
    res_scores = results_scores;
end

end

