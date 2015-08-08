function cluster_num_format(ax)
% CLUSTER_NUM_FORMAT format cluster number plot for use with SAVE_PNG
%
%   CLUSTER_NUM_FORMAT(ax)


for j = 1:length(ax.Children)
    if length(ax.Children(j).XData) == 1
        ax.Children(j).MarkerSize = 30;
        ax.Children(j).Color = [0 0 0];
    end
end
ax.XLabel = xlabel('Number Clusters');
ax.XLabel.FontSize = 16;
ax.YLabel.FontSize = 16;
ax.Title.FontSize = 18;
ax.FontSize = 14;
axis tight