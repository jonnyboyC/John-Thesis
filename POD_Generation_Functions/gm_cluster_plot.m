function gm_cluster_plot(h, modal_amp, gm_idx, gm_model, modes, k, direct, save_figures)
% Create 2D representation of clusters generated grom gaussian model

% Create new figure and axis handles
figure(h);
ax = newplot;

% Loop through and plot each cluster
hold(ax, 'on');
for i = 1:k
    figure(h);
    plot(ax, modal_amp(gm_idx == i, modes(1)), modal_amp(gm_idx == i, modes(2)), '.', ...
        'MarkerFaceColor', rand(1,3), 'MarkerEdgeColor', rand(1,3));
end

ezcontour(@(x1,x2)pdf(gm_model,[x1 x2]),[ax.XLim ax.YLim],100)
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
        'Clusters' filesep '2D_GM_Cluster'], save_figures{i});
end
end