function [scores_km, scores_gm] = classify_Gal(centers, gm_model, modal_amp, direct, km_stoch, gm_stoch)

make_plot = true;
save_figures = {};

h_km = figure;
h_gm = figure;

scores_km = zeros(length(modal_amp),1);
scores_gm = zeros(length(modal_amp),1);

for j = 1:length(modal_amp)
    gm_groups = cluster(gm_model, modal_amp{j,1}(:,2:end));
    km_groups = knnsearch(centers, modal_amp{j,1}(:,2:end));

    km_stoch_gal = gen_stochastic_matrix(h_km, 10, gm_groups, direct, make_plot, save_figures);
    gm_stoch_gal = gen_stochastic_matrix(h_gm, 10, km_groups, direct, make_plot, save_figures);
    
    scores_km(j) = norm(km_stoch_gal - km_stoch, 'fro');
    scores_gm(j) = norm(gm_stoch_gal - gm_stoch, 'fro');
end
scores_km = {scores_km};
scores_gm = {scores_gm};
end