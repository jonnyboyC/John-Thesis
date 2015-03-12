function handle = test_cluster_plot(modal_amp, cl_idx, cl_center, modes, k)
h = figure;
hold on;
for i = 1:k
    plot3(modal_amp(cl_idx == i, modes(1)), modal_amp(cl_idx == i, modes(2)), modal_amp(cl_idx == i, modes(3)), '.');
    plot3(cl_center(i,modes(1)), cl_center(i,modes(2)), cl_center(i,modes(3)), 'kx');
end
hold off;
if nargout == 1
    handle = h;
end
end