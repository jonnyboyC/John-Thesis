function h = gmm_cluster_plot(modal_amp, gmm, modes, k, direct, save_figures)
% GM_CLUSTER_PLOT generate a 2D for data clustered by a Gaussian mixture
% model
%
%   h = GM_CLUSTER_PLOT(modal_amp, gm, modes, k, direct, save_figures)

% Create new figure and axis handles
h = figure;
ax = newplot;
sub_ax = cell(4,1);

groups = gmm.groups;
models = gmm.models;

% Loop through and plot each cluster
hold(ax, 'on');
for i = 1:k
    figure(h);
    sub_ax{i} = plot(ax, modal_amp(groups == i, modes(1)), modal_amp(groups == i, modes(2)), '.', ...
        'MarkerFaceColor', rand(1,3), 'MarkerEdgeColor', rand(1,3));
end

% Plot prosterior probability contour lines
ezcontour(@(x1,x2)pdf(models,[x1 x2]),[ax.XLim, ax.YLim], 100);

for i = 1:k
    uistack(sub_ax{i}, 'top');
end

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