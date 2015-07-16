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
        sim_gm_groups = cluster(gm.models, integration{i}.modal_amp.(m{j}).(s{k})(:,2:end));
        sim_km_groups = knnsearch(km.centers, integration{i}.modal_amp.(m{j}).(s{k})(:,2:end));
        
        % Generate a stochastic matrix from transitions
        km_stoch_gal = gen_stochastic_matrix(num_clusters, sim_gm_groups);
        gm_stoch_gal = gen_stochastic_matrix(num_clusters, sim_km_groups);      
        
        % Determine the log probabiliy of observing a particular chain
        sim_prob_km = calc_probability(km.stoch, sim_km_groups);
        sim_prob_gm = calc_probability(gm.stoch, sim_gm_groups);
        
        frob_km.(m{j}).(s{k}) = norm(km_stoch_gal - km.stoch, 'fro');
        frob_gm.(m{j}).(s{k}) = norm(gm_stoch_gal - gm.stoch, 'fro');
        
        prob_km.(m{j}).(s{k}) = sim_prob_km - km.prob;
        prob_gm.(m{j}).(s{k}) = sim_prob_gm - gm.prob;
        
        plot_stochastic_matrix(km_stoch_gal, sim_km_groups, save_figures, direct, h_km);
        plot_stochastic_matrix(gm_stoch_gal, sim_gm_groups, save_figures, direct, h_gm);
        
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