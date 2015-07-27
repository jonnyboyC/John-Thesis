function mean_insta_fluc(results, name, image)

close all

X = results.X;
mean_W = results.mean_W;
flux_U = results.flux_U;
bnd_idx = results.bnd_idx;
bnd_X = results.bnd_X;
dimensions = results.dimensions;
uniform = results.uniform;

w = flow_comps(mean_W, X);
flux_UdX = derivatives(flux_U, bnd_idx, bnd_X, X, uniform, dimensions);
[x, u] = flow_comps_ip(X, flux_U);

if any(strcmp(u, 'w')) && any(strcmp(u, 'v'))
    flux_W.u = flux_UdX.w.y - flux_UdX.v.z;
end
if any(strcmp(u, 'u')) && any(strcmp(u, 'w'))
    flux_W.v = flux_UdX.u.z - flux_UdX.w.x;
end
if any(strcmp(u, 'v')) && any(strcmp(u, 'u'))
    flux_W.w = flux_UdX.v.x - flux_UdX.u.y;
end

% data.X = X;
for i = 1:2
    data.X.(x{i}) = X.(x{i});
end
data.bnd_idx = bnd_idx;
data.pod = reshape(mean_W.(w{1}), dimensions);
data.pod = data.pod;


h = figure;
plot_scalar_field(data);
ax = gca;
c = colorbar('southoutside');
c.Label.String  = 's^{-1}';
c.FontSize = 14;
ax.Title = title(['Mean Vorticity for ' name], 'FontSize', 18);
ax.XLabel = xlabel('X (mm)', 'FontSize', 16);
ax.YLabel = ylabel('Y (mm)', 'FontSize', 16);
saveas(h, ['D:\thesis\take2\mean' name], 'fig');

data.pod = reshape(mean_W.(w{1}), dimensions) + flux_W.(w{1})(:,:,image);
data.pod = data.pod;

h = figure;
plot_scalar_field(data);
ax = gca;
c = colorbar('southoutside');
c.Label.String = 's^{-1}';
c.FontSize = 14;
ax.Title = title(['Instantaneous Vorticity for ' name], 'FontSize', 18);
ax.XLabel = xlabel('X (mm)', 'FontSize', 16);
ax.YLabel = ylabel('Y (mm)', 'FontSize', 16);
saveas(h, ['D:\thesis\take2\insta' name], 'fig');


data.pod = flux_W.(w{1})(:,:,image);
data.pod = data.pod;


h = figure;
plot_scalar_field(data);
ax = gca;
c = colorbar('southoutside');
c.Label.String  = 's^{-1}';
c.FontSize = 14;
ax.Title = title(['Fluctuating Vorticity for ' name], 'FontSize', 18);
ax.XLabel = xlabel('X (mm)', 'FontSize', 16);
ax.YLabel = ylabel('Y (mm)', 'FontSize', 16);
saveas(h, ['D:\thesis\take2\fluc' name], 'fig');

