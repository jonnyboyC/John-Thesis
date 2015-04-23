function [pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac, l] = ...
    components(x, y, pod_u, pod_v, dimensions, vol_frac, num_modes, num_elem, bnd_idx, bnd_x, bnd_y)

% Old Methods will remove soon
% [xxi, yxi, xet, yet, aj] = metric2(x, y);    

% [pod_udx, pod_ud2x, pod_udy, pod_ud2y] = derivatives2(pod_u, x, y, dimensions);
% [pod_vdx, pod_vd2x, pod_vdy, pod_vd2y] = derivatives2(pod_v, x, y, dimensions);

% [a, b, c, d] = derivatives(pod_u, dimensions, ...
%     x, y, z, xxi, yxi, xet, yet, aj, bnd_idx);
% [e, f, g, h] = derivatives(pod_v, dimensions, ...
%     x, y, z, xxi, yxi, xet, yet, aj, bnd_idx);

% Calculate coefficients for for pod_u's & pod_v's derivatives
[pod_udx, pod_ud2x, pod_udy, pod_ud2y] = derivatives3(pod_u, bnd_idx, bnd_x, ...
                                            bnd_y, x, y, dimensions);
[pod_vdx, pod_vd2x, pod_vdy, pod_vd2y] = derivatives3(pod_v, bnd_idx, bnd_x, ...
                                            bnd_y, x, y, dimensions);

% Convert all matrix quantities to vectors
pod_udx     = reshape(pod_udx, num_elem, num_modes);
pod_udy     = reshape(pod_udy, num_elem, num_modes);
pod_ud2x    = reshape(pod_ud2x, num_elem, num_modes);
pod_ud2y    = reshape(pod_ud2y, num_elem, num_modes);
pod_vdx     = reshape(pod_vdx, num_elem, num_modes);
pod_vdy     = reshape(pod_vdy, num_elem, num_modes);
pod_vd2x    = reshape(pod_vd2x, num_elem, num_modes);
pod_vd2y    = reshape(pod_vd2y, num_elem, num_modes);

[pod_udx, pod_udy, pod_ud2x, pod_ud2y, pod_vdx, pod_vdy, pod_vd2x, pod_vd2y] = ...
    strip_boundaries(bnd_idx, pod_udx, pod_udy, pod_ud2x, pod_ud2y, pod_vdx, pod_vdy, pod_vd2x, pod_vd2y);

% Calculated Laplacian
d2pod_u = (pod_ud2x + pod_ud2y);
d2pod_v = (pod_vd2x + pod_vd2y);

[pod_u, pod_v, vol_frac] = strip_boundaries(bnd_idx, pod_u, pod_v, vol_frac);

% Find projected viscosity terms
cbu = inner_prod(d2pod_u, pod_u, vol_frac);
cbv = inner_prod(d2pod_v, pod_v, vol_frac);
l = cbu + cbv;

end