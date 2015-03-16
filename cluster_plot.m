function [handle, cax] = cluster_plot(ax, h, modal_amp, cl_idx, cl_center, modes, k)
if ax == 0
    h = figure;
    ax = newplot;
end
cla(ax);
hold(ax, 'on');
for i = 1:k
    plot3(ax, modal_amp(cl_idx == i, modes(1)), modal_amp(cl_idx == i, modes(2)), modal_amp(cl_idx == i, modes(3)), '.');
    plot3(ax, cl_center(i,modes(1)), cl_center(i,modes(2)), cl_center(i,modes(3)), 'kx');
end
hold(ax, 'off');
axis(ax, 'equal');
axis(ax, 'tight');
view(ax, [40, 60, 30]);
drawnow;

ax.XLabel.String = 'modal amplitude \phi_1';
ax.YLabel.String = 'modal amplitude \phi_2';
ax.ZLabel.String = 'modal amplitude \phi_3';
ax.Title.String = 'Clusters';

if nargout > 0
    handle = h;
end
if nargout > 1
    cax = ax;
end
end