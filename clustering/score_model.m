function [frob_km, frob_gm, like_km, like_gm, completed, km_steady, gmm_steady] = score_model(score_info)
%  SCORE_MODEL generate model scores based on method outline in John Chabot's
%  Thesis 2015 classify simulated data then apply score metrics. If
%  clusters were not precompiled from POD_Gen they will generate them here
%
%   [frob_km, frob_gm, like_km, like_gm, completed, km_steady, gmm_steady] 
%       = SCORE_MODEL(score_info)

% Unpack variables
km              = score_info.km;
gmm             = score_info.gmm;
integration     = score_info.integration;
tspan           = score_info.tspan;
modal_amp       = score_info.modal_amp;
num_clusters    = score_info.num_clusters;
int_time        = score_info.int_time;
num_cores       = score_info.num_cores;
modes           = score_info.modes;
custom          = score_info.custom;
outlier_mode    = score_info.outlier_mode;
direct          = score_info.direct;
multiplier      = score_info.multiplier;
MOD             = score_info.MOD;

if MOD
    modal_amp_til = score_info.modal_amp_til;
end

% Declaring loop variables
km_steady = struct;
gmm_steady = struct;
completed = struct;

% Free memory
clear score_info

% inputs for gen_stochoastic_matrix
save_figures = {};

% Get eddy visocity models
m = flow_comps(integration.modal_amp);
m_comps = flow_ncomps(integration.modal_amp);

h_km = figure;
h_gm = figure;

% Strip off mean mode from raw data
modal_amp = modal_amp(:,2:end);

% Select or generate clusters from data
if ~MOD
    
    % Post thesis we can change this back
    [gmm, km] = gen_clusters(modal_amp, modes, num_clusters, num_cores);
    
    %%%%%%%%%%%%%%%%%% TODO REVERSE AFTER THESIS %%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     if custom
%         [gmm, km] = gen_clusters(modal_amp, modes, num_clusters, num_cores);
%     else
%         if (length(modes) - 1) <= length(km)
%             idx = length(modes) - 1;
%             km = km{idx};
%             gmm = gmm{idx};
%         else
%             [gmm, km] = gen_clusters(modal_amp, modes, num_clusters, num_cores);
%         end
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
    

classify = true;

% Loop through each model calculating scores
for j = 1:m_comps
    
    s = flow_comps(integration.modal_amp.(m{j}));
    s_comps = flow_ncomps(integration.modal_amp.(m{j}));
    
    for k = 1:s_comps
        
        % If something happened the integration don't attempt to classify
        if isempty(integration.t.(m{j}).(s{k})) || isfield(integration.t.(m{j}).(s{k}), 'error')
            frob_km.(m{j}).(s{k}) = [];
            frob_gm.(m{j}).(s{k}) = [];
            
            like_km.(m{j}).(s{k}) = [];
            like_gm.(m{j}).(s{k}) = [];
            
            completed.(m{j}).(s{k}) = false;
            
            fprintf('%s %s took longer than %d seconds to integrate and was discarded\n', ...
                strrep((m{j}), '_', ' '), strrep((s{k}), '_', ' '), int_time);
            continue;
        end
        
        % Strip off first mode if not not modified
        if MOD
            modal_amp_model = integration.modal_amp.(m{j}).(s{k});
            [gmm, km] = gen_clusters(modal_amp_til.(m{j}).(s{k}), modes, ...
                num_clusters, num_cores, outlier_mode);
        else
            modal_amp_model = integration.modal_amp.(m{j}).(s{k})(:,2:end);
        end
        
        % Classify each point to a cluster
        [gmm_sim, km_sim] = classify_model(km, gmm, modal_amp_model, num_clusters, outlier_mode);
        
        % Generate a stochastic matrix from transitions
        km_sim = gen_stochastic_matrix(km_sim, num_clusters, multiplier, classify, outlier_mode);
        gmm_sim = gen_stochastic_matrix(gmm_sim, num_clusters, multiplier, classify, outlier_mode);
        
        % log probability of observing the simulate chain using MLE model
        km_sim.like = calc_likelihood(km.groups, km_sim.stoch);
        gmm_sim.like = calc_likelihood(gmm.groups, gmm_sim.stoch);
        
        % Frobenius norm of transition matrices
        frob_km.(m{j}).(s{k}) = norm(km_sim.stoch - km.stoch, 'fro');
        frob_gm.(m{j}).(s{k}) = norm(gmm_sim.stoch - gmm.stoch, 'fro');
        
        % Relative likelihood between observed chains
        like_km.(m{j}).(s{k}) = -(km_sim.like - km.like);
        like_gm.(m{j}).(s{k}) = -(gmm_sim.like- gmm.like);
        
        plot_stochastic_matrix(km_sim, save_figures, direct, h_km);
        plot_stochastic_matrix(gmm_sim, save_figures, direct, h_gm);
        
        km_steady.(m{j}).(s{k}) = km_sim.steady;
        gmm_steady.(m{j}).(s{k}) = gmm_sim.steady;
        
        % Temporary reporting
        if integration.t.(m{j}).(s{k})(end) == tspan(end)
            fprintf('%s %s finished got km_frob = %2.4f and gm_frob = %2.4f and km_prob = %2.4f and gm_prob = %2.4f\n', ...
                strrep((m{j}), '_', ' '), strrep((s{k}), '_', ' '),frob_km.(m{j}).(s{k}), ...
                frob_gm.(m{j}).(s{k}), like_km.(m{j}).(s{k}), like_gm.(m{j}).(s{k}));
            completed.(m{j}).(s{k}) = true;
        else
            fprintf('%s %s diverged got km_frob = %2.4f and gm_frob = %2.4f and km_prob = %2.4f and gm_prob = %2.4f\n', ...
                strrep((m{j}), '_', ' '), strrep((s{k}), '_', ' '),frob_km.(m{j}).(s{k}), ...
                frob_gm.(m{j}).(s{k}), like_km.(m{j}).(s{k}), like_gm.(m{j}).(s{k}));
            completed.(m{j}).(s{k}) = false;
        end
    end
end
end

