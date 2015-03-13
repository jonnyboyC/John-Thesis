function handle = cluster_plot(modal_amp, cl_idx, cl_center, modes, k)
h = figure;
ax = newplot;
hold(ax, 'on')
for i = 1:k
    plot3(ax, modal_amp(cl_idx == i, modes(1)), modal_amp(cl_idx == i, modes(2)), modal_amp(cl_idx == i, modes(3)), '.');
    plot3(ax, cl_center(i,modes(1)), cl_center(i,modes(2)), cl_center(i,modes(3)), 'kx');
end
hold(ax, 'off');
axis(ax, 'equal');
axis(ax, 'tight');
view(ax, [40, 60, 30])


[vertex, center] = voronoin(cl_center);


if nargout == 1
    handle = h;
end
end