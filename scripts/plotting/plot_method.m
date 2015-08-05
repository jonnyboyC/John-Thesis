function plot_method(bnd_idx, bnd_X, symbols)

h = figure;
ax = gca;

ax.XLim = [1 size(bnd_idx,1) + 1];
ax.XTick = 1:size(bnd_idx,1) + 1;
ax.YLim = [1 size(bnd_idx,2) + 1];
ax.YTick = 1:size(bnd_idx,2) + 1;
hold(ax, 'on');
grid(ax, 'on');

for i = 1:size(bnd_idx,1)
    for j = 1:size(bnd_idx,2)
        if bnd_idx(i,j) == -1
            patch = fill([i, i, i+1, i+1], [j, j+1, j+1, j], 'k');
            patch.FaceAlpha = 0.3;
        end
        if bnd_idx(i,j) == 0 && bnd_X.x(i,j) == 0 && bnd_X.y(i,j) == 0
            patch = fill([i, i, i+1, i+1], [j, j+1, j+1, j], 'r');
            patch.FaceAlpha = 0.15;
        end
        if bnd_idx(i,j) == 0 && (~bnd_X.x(i,j) == 0 || ~bnd_X.y(i,j) == 0)
            patch = fill([i, i, i+1, i+1], [j, j+1, j+1, j], 'k');
            patch.FaceAlpha = 0.15;
        end
        text(i+0.1,j+0.5, symbols{i,j});
    end
end