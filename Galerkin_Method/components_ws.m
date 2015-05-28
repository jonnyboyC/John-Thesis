function [pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac] = ...
        components_ws(x, y, pod_u, pod_v, dimensions, vol_frac, num_modes, bnd_idx, bnd_x, bnd_y)
% COMPONENT_WS, calculate the derivative terms needed to calculate the
% galerkin coefficients of the weak formulation
% 
% type edit COMPONENT for input details
    

% Calculate coefficients for for pod_u's & pod_v's derivatives
[pod_udx, pod_udy] = derivatives(pod_u, bnd_idx, bnd_x, bnd_y, x, y, dimensions);
[pod_vdx, pod_vdy] = derivatives(pod_v, bnd_idx, bnd_x, bnd_y, x, y, dimensions);

pod_udx = reshape(pod_udx, [], num_modes);
pod_udy = reshape(pod_udy, [], num_modes);
pod_vdx = reshape(pod_vdx, [], num_modes);
pod_vdy = reshape(pod_vdy, [], num_modes);

[pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac] = ...
    strip_boundaries(bnd_idx, pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac);

end