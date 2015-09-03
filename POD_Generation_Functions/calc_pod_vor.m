function [pod_W, mean_W] = calc_pod_vor(U, mean_U, dimensions, ...
                                    bnd_idx, bnd_X, X)
% CALC_POD_VOR determine the vorticity component of each POD mode
%
% [pod_vor, mean_vor] = CALC_POD_VOR(u, v, mean_u, mean_v, dimensions,
% cutoff, bnd_idx, bnd_x, bnd_, x, y)

% Calculate vorticity values using galerkin derivatives function
[~, u] = flow_comps_ip(X, U);

% Calculate dervatives
UdX = derivatives(U, bnd_idx, bnd_X, X, dimensions); %
mean_UdX = derivatives(mean_U, bnd_idx, bnd_X, X, dimensions);


% Calculate each vorticity component based on the present flow components
if any(strcmp(u, 'w')) && any(strcmp(u, 'v'))
    pod_W.u = UdX.w.y - UdX.v.z;
    mean_W.u = mean_UdX.w.y - mean_UdX.v.z;
end
if any(strcmp(u, 'u')) && any(strcmp(u, 'w'))
    pod_W.v = UdX.u.z - UdX.w.x;
    mean_W.v = mean_UdX.u.z - mean_UdX.w.x;
end
if any(strcmp(u, 'v')) && any(strcmp(u, 'u'))
    pod_W.w = UdX.v.x - UdX.u.y;
    mean_W.w = mean_UdX.v.x - mean_UdX.u.y;
end
end