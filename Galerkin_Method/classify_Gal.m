function [frob_km, frob_gm, prob_km, prob_gm] = classify_Gal(km, gm, integration, i, direct)
 
% inputs for gen_stochoastic_matrix
make_plot = false;
save_figures = {};

% Get eddy visocity models
e = flow_comps(integration{i}.modal_amp);
e_comps = flow_ncomps(integration{i}.modal_amp);

h_km = figure;
h_gm = figure;

% Loop through each model calculating scores
for j = 1:e_comps
    
    s = flow_comps(integration{i}.modal_amp.(e{j}));
    s_comps = flow_ncomps(integration{i}.modal_amp.(e{j}));
    
    for k = 1:s_comps
        sim_gm_groups = cluster(gm.models, integration{i}.modal_amp.(e{j}).(s{k})(:,2:end));
        sim_km_groups = knnsearch(km.centers, integration{i}.modal_amp.(e{j}).(s{k})(:,2:end));
        
        km_stoch_gal = gen_stochastic_matrix(h_km, 10, sim_gm_groups, direct, make_plot, save_figures);
        gm_stoch_gal = gen_stochastic_matrix(h_gm, 10, sim_km_groups, direct, make_plot, save_figures);
        
        frob_km.(e{j}).(s{k}) = norm(km_stoch_gal - km.stoch, 'fro');
        frob_gm.(e{j}).(s{k}) = norm(gm_stoch_gal - gm.stoch, 'fro');
        
        sim_prob_km = calc_probability(km.stoch, sim_km_groups);
        sim_prob_gm = calc_probability(gm.stoch, sim_gm_groups);
        
        prob_km.(e{j}).(s{k}) = exp(sim_prob_km - km.prob);
        prob_gm.(e{j}).(s{k}) = exp(sim_prob_gm - gm.prob);
    end
end
end