function [pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac] = ...
        components_ws(x, y, pod_u, pod_v, dimensions, vol_frac, num_modes, bnd_idx, bnd_x, bnd_y)

% TODO figure out what is really being calculated here
% [xxi, yxi, xet, yet, aj] = metric2(x, y);

% Calculate coefficients for for pod_u's & pod_v's derivatives
% [pod_udx, ~, pod_udy, ~] = derivatives(pod_u, dimensions, ...
%     z, xxi, yxi, xet, yet, aj, bnd_idx);
% [pod_vdx, ~, pod_vdy, ~] = derivatives(pod_v, dimensions, ...
%     z, xxi, yxi, xet, yet, aj, bnd_idx);

% [pod_udx, ~, pod_udy, ~] = derivatives2(pod_u, x, y, dimensions);
% [pod_vdx, ~, pod_vdy, ~] = derivatives2(pod_v, x, y, dimensions);

[pod_udx, ~, pod_udy, ~] = derivatives3(pod_u, bnd_idx, bnd_x, bnd_y, x, y, dimensions);
[pod_vdx, ~, pod_vdy, ~] = derivatives3(pod_v, bnd_idx, bnd_x, bnd_y, x, y, dimensions);

pod_udx = reshape(pod_udx, [], num_modes);
pod_udy = reshape(pod_udy, [], num_modes);
pod_vdx = reshape(pod_vdx, [], num_modes);
pod_vdy = reshape(pod_vdy, [], num_modes);

[pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac] = ...
    strip_boundaries(bnd_idx, pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac);

end