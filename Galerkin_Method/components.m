function [pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac, l] = ...
    components(x, y, mean_u, mean_v, pod_u, pod_v, dimensions, vol_frac, num_modes, num_elem, bnd_idx, z)

if num_modes > 30
    if isempty(gcp);
        parpool('local', 3);
    end
end

% TODO figure out what is really being calculated here
[xxi, yxi, xet, yet, aj] = metric2(x, y);

pod_u = cat(3,mean_u, pod_u);
pod_v = cat(3,mean_v, pod_v);

% Calculate coefficients for for pod_u's & pod_v's derivatives
[pod_udx, pod_ud2x, pod_udy, pod_ud2y] = derivatives(pod_u, dimensions, ...
    z, xxi, yxi, xet, yet, aj, bnd_idx);
[pod_vdx, pod_vd2x, pod_vdy, pod_vd2y] = derivatives(pod_v, dimensions, ...
    z, xxi, yxi, xet, yet, aj, bnd_idx);

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

if ~isempty(gcp)
    delete(gcp);
end
end