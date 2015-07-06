function [km, gm] = cluster_POD(modal_amp, num_clusters, direct, save_figures)
% CLUSTER_POD return clusters of the POD modal amplitudes by k-mean and
% gaussian mixture model type edit to find list of imputs and outputs

% Flag
make_plot = true;

% k-means replicate
k_replicate = 50;
gm_replicate = 10;

% Cluster dimensions of clusters
cluster_range = 2:40;
cluster_modes = [1,2];

% Prefil groups and centers
gm = cell(length(cluster_range),1);
km = cell(length(cluster_range),1);

% Clustering options
options = statset('UseParallel', 1);
gm_options = statset('MaxIter', 1000);

% ready figures
h_km = figure;
h_gm = figure;
h_stoch = figure;
h_stoch2 = figure;

% Calculate cluster centers and stochastic matrices
for i = 1:length(cluster_range);
    % Cluster based on k-mean
    [km{i}.groups, km{i}.centers] = kmeans(modal_amp(:,1:cluster_range(i)), ...
        num_clusters, 'Replicates', k_replicate, 'Options', options);
    
    % Create gaussian mixture model
    gm{i}.models = fitgmdist(modal_amp(:,1:cluster_range(i)), num_clusters, ...
        'Replicates', gm_replicate, 'SharedCov', true, 'Options', gm_options);
    
    % Cluster based on gaussian mixture model
    gm{i}.groups = cluster(gm{i}.models, modal_amp(:,1:cluster_range(i)));
    
    % in 2D phase space plot clusters
    if i == 1 && ~isempty(save_figures)
        cluster_plot(h_km, modal_amp, km{i}.groups, km{i}.centers, cluster_modes, num_clusters, ...
            direct, save_figures);
        gm_cluster_plot(h_gm, modal_amp, gm{i}.groups, gm{i}.models, cluster_modes, ...
            num_clusters, direct, save_figures);
    end
    
    % Generate stochastic matrices for both clustering methods
    km{i}.stoch = gen_stochastic_matrix(h_stoch, num_clusters, km{i}.groups, direct, make_plot, save_figures);
    km{i}.prob = calc_probability(km{i}.stoch, km{i}.groups);
    gm{i}.stoch = gen_stochastic_matrix(h_stoch2, num_clusters, gm{i}.groups, direct, make_plot, save_figures);
    gm{i}.prob = calc_probability(gm{i}.stoch, gm{i}.groups);
end

end