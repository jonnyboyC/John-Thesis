function [km_stoch, gm_stoch, gm_models, gm_groups, km_groups, centers] = ...
    cluster_POD(modal_amp, num_clusters, direct, save_figures)

% Flag
make_plot = true;

% Cluster dimensions of clusters
cluster_range = 2:40;
cluster_modes = [1,2];

% Prefil groups and centers
centers = cell(length(cluster_range),1);
gm_models = cell(length(cluster_range),1);
gm_groups = cell(length(cluster_range),1);
km_stoch = cell(length(cluster_range),1);
gm_stoch = cell(length(cluster_range),1);

% Clustering options
options = statset('UseParallel', 1);
gm_options = statset('MaxIter', 1000);

% ready figures
h_km = figure;
h_gm    = figure;
h_stoch = figure;
h_stoch2= figure;

% Calculate cluster centers and stochastic matrices
for i = 1:length(cluster_range);
    % Cluster based on k-mean
    [km_groups, centers{i}] = kmeans(modal_amp(:,1:cluster_range(i)), ...
        num_clusters, 'Replicates', 10, 'Options', options);
    
    % Create gaussian mixture model
    gm_models{i} = fitgmdist(modal_amp(:,1:cluster_range(i)), num_clusters, ...
        'Replicates', 10, 'SharedCov', true, 'Options', gm_options);
    
    % Cluster based on gaussian mixture model
    gm_groups{i} = cluster(gm_models{i}, modal_amp(:,1:cluster_range(i)));
    
    % in 2D phase space plot clusters
    if i == 1 && ~isempty(save_figures)
        cluster_plot(h_km, modal_amp, km_groups, centers{i}, cluster_modes, num_clusters, ...
            direct, save_figures);
        gm_cluster_plot(h_gm, modal_amp, gm_groups{i}, gm_models{i}, cluster_modes, ...
            num_clusters, direct, save_figures);
    end
    
    % Generate stochastic matrices for both clustering methods
    km_stoch{i} = gen_stochastic_matrix(h_stoch, km_groups, direct, make_plot, save_figures);
    gm_stoch{i} = gen_stochastic_matrix(h_stoch2, gm_groups{i}, direct, make_plot, save_figures);
end

end