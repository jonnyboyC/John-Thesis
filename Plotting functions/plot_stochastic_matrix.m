function h = plot_stochastic_matrix(model, save_figures, direct, h)
% PLOT_STOCHASTIC_MATRIX plot a predicted stochastic matrix for a clustered
% model
%
%   h = PLOT_STOCHASTIC_MATRIX(model, save_figures, direct, h)

if nargin == 3
    h = figure;
end

if ~isvalid(h)
    h = figure;
else
    figure(h);
end

groups = model.groups;
stoch = model.stoch;

clusters = length(stoch);

with_outlier = zeros(clusters+1);
with_outlier(1:clusters,1:clusters) = stoch;

% plot transition matrix
pcolor(with_outlier);
ax = gca;
colorbar;

ax.Title = title('stochastic matrix');
ax.XLabel = xlabel('Next State');
ax.YLabel = ylabel('Current State');

Label = cell(clusters, 1);
for i = 1:clusters-1
    Label{i} = num2str(i);
end
Label{end} = 'Out';

ax.XTick = (1:clusters) + 0.5;
ax.YTick = (1:clusters) + 0.5;

ax.XTickLabel = Label;
ax.YTickLabel = Label;

hold(ax, 'on');
for i = 1:clusters
    for j = 1:clusters
        prob = stoch(i,j);
        rectangle('Position', [j+0.5-prob/2, i+0.5-prob/2, prob, prob], ...
            'Curvature', [1,1]);
    end
end
hold(ax, 'off');

if ~isempty(save_figures)
    for i = 1:length(save_figures)
        saveas(h, [direct filesep 'Figures' filesep 'POD' filesep ...
            'Clusters' filesep 'Stochastic_matirx_modes', num2str(size(groups,2))], save_figures{i});
    end
end

drawnow;