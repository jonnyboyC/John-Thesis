function h = km_cluster_plot(modal_amp, km, cluster_modes, k, direct, save_figures)
% Create 2D representation of clusters centers

% Create new figure and axis handles
h = figure;
ax = newplot;

centers = km.centers;
groups = km.groups;

% Loop through and plot each cluster
hold(ax, 'on');
for i = 1:k
    figure(h);
    plot(ax, modal_amp(groups == i, cluster_modes(1)), modal_amp(groups == i, cluster_modes(2)), '.', ...
        'MarkerFaceColor', rand(1,3), 'MarkerEdgeColor', rand(1,3));
    plot(ax, centers(i,cluster_modes(1)), centers(i,cluster_modes(2)), 'kx');
end

% Generate voronoi cells
voronoi(ax, centers(1:k,cluster_modes(1)), centers(1:k,cluster_modes(2)))

% Formatting
hold(ax, 'off');
axis(ax, 'equal');
axis(ax, 'tight');
drawnow;

% Add axis labels
ax.XLabel.String = 'modal amplitude \phi_1';
ax.YLabel.String = 'modal amplitude \phi_2';
ax.Title.String = 'Clusters';

for i = 1:length(save_figures)
    saveas(h, [direct filesep 'Figures' filesep 'POD' filesep ...
        'Clusters' filesep '2D_Cluster'], save_figures{i});
end
end