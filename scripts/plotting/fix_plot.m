files = ls;
files = files(3:end, :);
for i = 1:size(files,1)
    h = hgload(files(i,:));
    ax = gca;
    ax.CLim = [-max(abs(ax.CLim)), max(abs(ax.CLim))];
    saveas(h, files(i,:), 'fig')
end