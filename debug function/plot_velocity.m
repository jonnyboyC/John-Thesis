function plot_velocity(u, v, x, y, i, dimensions, mean_u, mean_v)
if nargin == 8
    u = u(:,i) - mean_u;
    v = v(:,i) - mean_v;
else
    u = u(:,i);
    v = v(:,i);
end
mag = reshape(sqrt(u.^2 + v.^2), dimensions);
pcolor(x, y, mag)
shading interp
hold on
quiver(x, y, reshape(u, dimensions), reshape(v, dimensions), 'k');
hold off
pause(.2);
end