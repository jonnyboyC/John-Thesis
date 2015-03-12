function [pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac] = ...
        components_ws(x, y, pod_u, pod_v, dimensions, num_modes, vol_frac, bnd_idx)

if num_modes > 30
    if isempty(gcp);
        parpool('local', 3);
    end
end

% TODO figure out what is really being calculated here
[xxi, yxi, xet, yet, aj] = metric2(x, y);

% Calculate coefficients for for pod_u's & pod_v's derivatives
[pod_udx, ~, pod_udy, ~] = derivatives(pod_u, dimensions, ...
    z, xxi, yxi, xet, yet, aj, bnd_idx);
[pod_vdx, ~, pod_vdy, ~] = derivatives(pod_v, dimensions, ...
    z, xxi, yxi, xet, yet, aj, bnd_idx);

[pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac] = ...
    strip_boundaries(bnd_idx, pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac);
end