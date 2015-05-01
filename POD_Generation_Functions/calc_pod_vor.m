function pod_vor = calc_pod_vor(u, v, dimensions, cutoff, bnd_idx, bnd_x, bnd_y, x, y)
% Calculate vorticity values using galerkin derivatives function
pod_vor = zeros(dimensions(1), dimensions(2), cutoff);

% Calculate dervatives
[~, udy] = derivatives(u, bnd_idx, bnd_x, bnd_y, x, y, dimensions);
[vdx, ~] = derivatives(v, bnd_idx, bnd_x, bnd_y, x, y, dimensions);

for i = 1:cutoff
    pod_vor(:,:,i) = vdx(:,:,i)-udy(:,:,i);
end
end