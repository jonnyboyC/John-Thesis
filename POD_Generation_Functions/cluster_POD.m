function [km, gm] = cluster_POD(modal_amp, num_clusters, num_cores, direct, save_figures)
% CLUSTER_POD return clusters of the POD modal amplitudes by k-mean and
% gaussian mixture model type edit to find list of imputs and outputs

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
if num_cores ~= 1
    options = statset('UseParallel', true);
else
    options = statset('UseParallel', false);
end
gm_options = statset('MaxIter', 1000);

% ready figures
h_stoch = figure;
h_stoch2 = figure;

valid = false;
multiplier = 1;

% Calculate cluster centers and stochastic matrices
for i = 1:length(cluster_range);
    % Cluster based on k-mean
    [km{i}.groups, km{i}.centers] = kmeans(modal_amp(:,1:cluster_range(i)), ...
        num_clusters, 'Replicates', k_replicate, 'Options', options);
    
    % Create gaussian attempt to use independent covariance matrix first
    try 
        gm{i}.models = fitgmdist(modal_amp(:,1:cluster_range(i)), num_clusters, ...
            'Replicates', gm_replicate, 'SharedCov', false, 'Options', gm_options);
    catch error
        if (strcmp(error.identifer, 'stats:gmdistribution:IllCondCov')
            gm{i}.models = fitgmdist(modal_amp(:,1:cluster_range(i)), num_clusters, ...
                'Replicates', gm_replicate, 'SharedCov', true, 'Options', gm_options);
        else
            rethrow(error);
        end
    end
    
    % Cluster based on gaussian mixture model
    gm{i}.groups = cluster(gm{i}.models, modal_amp(:,1:cluster_range(i)));
    
    % in 2D phase space plot clusters
    if i == 1 && ~isempty(save_figures)
        cluster_plot(modal_amp, km{i}.groups, km{i}.centers, cluster_modes, num_clusters, ...
            direct, save_figures);
        gm_cluster_plot(modal_amp, gm{i}.groups, gm{i}.models, cluster_modes, ...
            num_clusters, direct, save_figures);
    end
    
    % Generate stochastic matrices for both clustering methods
    km{i}.stoch = gen_stochastic_matrix(num_clusters, km{i}.groups, multiplier, valid);
    gm{i}.stoch = gen_stochastic_matrix(num_clusters, gm{i}.groups, multiplier, valid);
    km{i}.stat = stationary(km{i}.stoch);
    gm{i}.stat = stationary(gm{i}.stoch);
    
    gm{i}.prob = calc_probability(gm{i}.stoch, gm{i}.groups);
    km{i}.prob = calc_probability(km{i}.stoch, km{i}.groups);
    
    plot_stochastic_matrix(km{i}.stoch, km{i}.groups, save_figures, direct, h_stoch);
    plot_stochastic_matrix(gm{i}.stoch, gm{i}.groups, save_figures, direct, h_stoch2);

end

end