function [udx, udy, vdx, vdy, pod_udx, pod_udy, pod_vdx, pod_vdy] = ...
        components_ws(x, y, mean_u, mean_v, pod_u, pod_v, dimensions, num_modes)

if num_modes > 30
    if isempty(gcp);
        parpool('local', 3);
    end
end

% TODO figure out what is really being calculated here
[xxi, yxi, xet, yet, aj] = metric2(x, y);
    
% Calculate coefficients for u's & v's derivatives
[udx, ~, udy, ~] = derivatives(mean_u, dimensions, z, xxi, yxi,...
    xet, yet, aj, bnd_idx);
[vdx, ~, vdy, ~] = derivatives(mean_v, dimensions, z, xxi, yxi,...
    xet, yet, aj, bnd_idx);

% Calculate coefficients for for pod_u's & pod_v's derivatives
[pod_udx, ~, pod_udy, ~] = derivatives(pod_u, dimensions, ...
    z, xxi, yxi, xet, yet, aj, bnd_idx);
[pod_vdx, ~, pod_vdy, ~] = derivatives(pod_v, dimensions, ...
    z, xxi, yxi, xet, yet, aj, bnd_idx);

end