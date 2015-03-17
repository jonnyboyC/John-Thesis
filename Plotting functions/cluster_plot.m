function cluster_plot(modal_amp, cl_idx, cl_center, modes, k, direct, save_figures)
% Create 2D representation of clusters centers

% Create new figure and axis handles
h = figure;
ax = newplot;

% Loop through and plot each cluster
hold(ax, 'on');
for i = 1:k
    figure(h);
    plot(ax, modal_amp(cl_idx == i, modes(1)), modal_amp(cl_idx == i, modes(2)), '.', ...
        'MarkerFaceColor', rand(1,3), 'MarkerEdgeColor', rand(1,3));
    plot(ax, cl_center(i,modes(1)), cl_center(i,modes(2)), 'kx');
end

% Generate voronoi cells
voronoi(ax, cl_center(1:k,modes(1)), cl_center(1:k,modes(2)))

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