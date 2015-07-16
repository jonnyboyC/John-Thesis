function [frob_km, frob_gm, prob_km, prob_gm, completed] = classify_Gal(km, gm, integration, tspan, num_clusters, int_time, i, direct)
 
% inputs for gen_stochoastic_matrix
save_figures = {};

% Get eddy visocity models
m = flow_comps(integration{i}.modal_amp);
m_comps = flow_ncomps(integration{i}.modal_amp);

h_km = figure;
h_gm = figure;

% Loop through each model calculating scores
for j = 1:m_comps
    
    s = flow_comps(integration{i}.modal_amp.(m{j}));
    s_comps = flow_ncomps(integration{i}.modal_amp.(m{j}));
    
    for k = 1:s_comps
        
        % If something happened the integration don't attempt to classify
        if isempty(integration{i}.t.(m{j}).(s{k})) || isfield(integration{i}.t.(m{j}).(s{k}), 'error')
            frob_km.(m{j}).(s{k}) = [];
            frob_gm.(m{j}).(s{k}) = [];
            
            prob_km.(m{j}).(s{k}) = [];
            prob_gm.(m{j}).(s{k}) = [];
            
            fprintf('%s %s took longer than %d seconds to integrate and was discarded\n', ...
                strrep((m{j}), '_', ' '), strrep((s{k}), '_', ' '), int_time);
            continue;
        end
        
        % Classify each point to a cluster
        sim_km_groups = knnsearch(km.centers, integration{i}.modal_amp.(m{j}).(s{k})(:,2:end));
        sim_gm_groups = cluster(gm.models, integration{i}.modal_amp.(m{j}).(s{k})(:,2:end));
        
        % Generate a stochastic matrix from transitions
        km_MLE = gen_stochastic_matrix(num_clusters, sim_km_groups);
        gm_MLE = gen_stochastic_matrix(num_clusters, sim_gm_groups);
        
        % log probability of observing the simulated chain using emp model
        km_prob_ALE = calc_probability(km.stoch, sim_km_groups);
        gm_prob_ALE = calc_probability(gm.stoch, sim_gm_groups);
        
        % log probability of observing the simulate chain using MLE model
        km_prob_MLE = calc_probability(km_MLE, sim_km_groups);
        gm_prob_MLE = calc_probability(gm_MLE, sim_gm_groups);
        
        % Frobenius norm of transition matrices
        frob_km.(m{j}).(s{k}) = norm(km_MLE - km.stoch, 'fro');
        frob_gm.(m{j}).(s{k}) = norm(gm_MLE - gm.stoch, 'fro');
        
        % Relative likelihood between observed chains
        prob_km.(m{j}).(s{k}) = km_prob_ALE - km_prob_MLE;
        prob_gm.(m{j}).(s{k}) = gm_prob_ALE - gm_prob_MLE;
        
        plot_stochastic_matrix(km_MLE, sim_km_groups, save_figures, direct, h_km);
        plot_stochastic_matrix(gm_MLE, sim_gm_groups, save_figures, direct, h_gm);
        
        % Temporary reporting
        if integration{i}.t.(m{j}).(s{k})(end) == tspan(end)
            fprintf('%s %s finished got km_frob = %2.4f and gm_frob = %2.4f and km_prob = %2.4f and gm_prob = %2.4f\n', ...
                strrep((m{j}), '_', ' '), strrep((s{k}), '_', ' '),frob_km.(m{j}).(s{k}), ...
                frob_gm.(m{j}).(s{k}), prob_km.(m{j}).(s{k}), prob_gm.(m{j}).(s{k}));
            completed.(m{j}).(s{k}) = true;
        else
            fprintf('%s %s diverged got km_frob = %2.4f and gm_frob = %2.4f and km_prob = %2.4f and gm_prob = %2.4f\n', ...
                strrep((m{j}), '_', ' '), strrep((s{k}), '_', ' '),frob_km.(m{j}).(s{k}), ...
                frob_gm.(m{j}).(s{k}), prob_km.(m{j}).(s{k}), prob_gm.(m{j}).(s{k}));
            completed.(m{j}).(s{k}) = false;
        end
    end
end
end