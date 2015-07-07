function h = plot_stochastic_matrix(stoc_matrix, groups, save_figures, direct, h)
% Plot a given stochastic matrix

if nargin == 4
    h = figure;
end

if ~isvalid(h)
    h = figure;
end

% plot transition matrix
pcolor(stoc_matrix);
ax = gca;
colorbar;
ax.Title = title('stochastic matrix');
ax.XLabel = xlabel('Next State');
ax.YLabel = ylabel('Current State');

if ~isempty(save_figures)
    for i = 1:length(save_figures)
        saveas(h, [direct filesep 'Figures' filesep 'POD' filesep ...
            'Clusters' filesep 'Stochastic_matirx_modes', num2str(size(groups,2))], save_figures{i});
    end
end

drawnow;