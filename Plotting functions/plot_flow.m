function ax = plot_flow(x, y, u, v, bnd_idx, image)

if size(size(u)', 1) == 2 && size(size(v)',1) == 2
    magnitude = sqrt(u.^2 + v.^2);
else
    magnitude = sqrt(u(:,:,image).^2 + v(:,:,image).^2);
end
data.pod = magnitude;
data.x = x;
data.y = y;
handle = figure;

[~, ax] = Plottec2(data, handle, bnd_idx);
ax.Title.String = ['Instanteous Flow Visualisation for Snapshot ' num2str(image)];
ax.XLabel.String = 'x/D';
ax.YLabel.String = 'y/D';
hold on

spacing_x = ceil(size(x, 1)/50);
spacing_y = ceil(size(x, 2)/50);

short_x = x(1:spacing_x:end, 1:spacing_y:end);
short_y = y(1:spacing_x:end, 1:spacing_y:end);
short_u = u(1:spacing_x:end, 1:spacing_y:end, image);
short_v = v(1:spacing_x:end, 1:spacing_y:end, image);

quiver(short_x, short_y, short_u, short_v, 'color', [0 0 0]);
hold off


if nargout
    h = handle;
end