function [pod_vor, mean_vor] = calc_pod_vor(u, v, mean_u, mean_v, dimensions, ...
    cutoff, bnd_idx, bnd_x, bnd_y, x, y)
% CALC_POD_VOR determine the vorticity component of each POD mode
%
% [pod_vor, mean_vor] = CALC_POD_VOR(u, v, mean_u, mean_v, dimensions,
% cutoff, bnd_idx, bnd_x, bnd_, x, y)

% Calculate vorticity values using galerkin derivatives function
pod_vor = zeros(dimensions(1), dimensions(2), cutoff);

% Calculate dervatives
[~, udy] = derivatives(u, bnd_idx, bnd_x, bnd_y, x, y, dimensions);
[vdx, ~] = derivatives(v, bnd_idx, bnd_x, bnd_y, x, y, dimensions);

for i = 1:cutoff
    pod_vor(:,:,i) = vdx(:,:,i)-udy(:,:,i);
end

[~, udy] = derivatives(mean_u, bnd_idx, bnd_x, bnd_y, x, y, dimensions);
[vdx, ~] = derivatives(mean_v, bnd_idx, bnd_x, bnd_y, x, y, dimensions);

mean_vor = vdx(:,:,1) - udy(:,:,1);
end