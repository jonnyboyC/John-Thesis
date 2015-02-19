function [udx, udy, vdx, vdy, pod_udx, pod_udy, pod_vdx, pod_vdy, mean_u, mean_v, pod_u, pod_v, vol_frac, l_dot, l] = ...
    components(x, y, mean_u, mean_v, pod_u, pod_v, dimensions, vol_frac, num_modes, num_elem, bnd_idx, z)

if num_modes > 30
    if isempty(gcp);
        parpool('local', 3);
    end
end

% TODO figure out what is really being calculated here
[xxi, yxi, xet, yet, aj] = metric2(x, y);

% Calculate coefficients for u's & v's derivatives
[udx, ud2x, udy, ud2y] = derivatives(mean_u, dimensions, z, xxi, yxi,...
    xet, yet, aj, bnd_idx);
[vdx, vd2x, vdy, vd2y] = derivatives(mean_v, dimensions, z, xxi, yxi,...
    xet, yet, aj, bnd_idx);

[mean_u, mean_v] = strip_boundaries(bnd_idx, mean_u, mean_v);

% Calculate coefficients for for pod_u's & pod_v's derivatives
[pod_udx, pod_ud2x, pod_udy, pod_ud2y] = derivatives(pod_u, dimensions, ...
    z, xxi, yxi, xet, yet, aj, bnd_idx);
[pod_vdx, pod_vd2x, pod_vdy, pod_vd2y] = derivatives(pod_v, dimensions, ...
    z, xxi, yxi, xet, yet, aj, bnd_idx);

% Convert all matrix quantities to vectors
udx     = reshape(udx, num_elem, 1);
udy     = reshape(udy, num_elem, 1);
ud2x    = reshape(ud2x, num_elem, 1);
ud2y    = reshape(ud2y, num_elem, 1);
vdx     = reshape(vdx, num_elem, 1);
vdy     = reshape(vdy, num_elem, 1);
vd2x    = reshape(vd2x, num_elem, 1);
vd2y    = reshape(vd2y, num_elem, 1);

[udx, udy, ud2x, ud2y, vdx, vdy, vd2x, vd2y] = ...
    strip_boundaries(bnd_idx, udx, udy, ud2x, ud2y, vdx, vdy, vd2x, vd2y);

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

% Need to check this out, should it be distance formula?
d2u = (ud2x + ud2y);
d2v = (vd2x + vd2y);
d2pod_u = (pod_ud2x + pod_ud2y);
d2pod_v = (pod_vd2x + pod_vd2y);

[pod_u, pod_v, vol_frac] = strip_boundaries(bnd_idx, pod_u, pod_v, vol_frac);

% Find inner product between 2nd for mean velocity and pod modes
clu = inner_prod(d2u, pod_u, vol_frac);
clv = inner_prod(d2v, pod_v, vol_frac);
l_dot = clu + clv;

% Find inner product between 2nd derivative of pod modes and pod modes
cbu = inner_prod(d2pod_u, pod_u, vol_frac);
cbv = inner_prod(d2pod_v, pod_v, vol_frac);
l = cbu + cbv;

if ~isempty(gcp)
    delete(gcp);
end
end