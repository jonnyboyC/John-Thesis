function [km, gmm] = cluster_POD(modal_amp, num_clusters, num_cores, garbage_mode, direct, save_figures)
% CLUSTER_POD return clusters of the POD modal amplitudes by k-mean and
% gaussian mixture. Generates modes for 2-40 modes, with plots
%
% [km, gm] = CLUSTER_POD(modal_amp, num_clusters, num_cores, direct, ...
%   save_figures)

% Cluster dimensions of clusters
cluster_range = 2:40;
cluster_modes = [1,2];

% Prefil groups and centers
gmm = cell(length(cluster_range),1);
km = cell(length(cluster_range),1);

% ready figures
h_stoch = figure;
h_stoch2 = figure;

% Calculate cluster centers and stochastic matrices
for i = 1:length(cluster_range);

    % Modes to be clustered
    modes = 1:cluster_range(i);
    
    % Generate clusters for the empirical data
    [gmm{i}, km{i}] = gen_clusters(modal_amp, modes, num_clusters, num_cores, garbage_mode);
    
    % in 2D phase space plot clusters
    if i == 1 && ~isempty(save_figures)
        cluster_plot(modal_amp, km{i}, cluster_modes, num_clusters, ...
            direct, save_figures);
        gm_cluster_plot(modal_amp, gmm{i}, cluster_modes, num_clusters, ...
            direct, save_figures);
    end
    
    % Get stationary distrubtion
    km{i} = stationary(km{i});
    gmm{i} = stationary(gmm{i});
    
    % Plot transition matrices
    plot_stochastic_matrix(km{i}, save_figures, direct, h_stoch);
    plot_stochastic_matrix(gmm{i}, save_figures, direct, h_stoch2);
end

end