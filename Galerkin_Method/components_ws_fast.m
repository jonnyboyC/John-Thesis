function [pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac] = ...
        components_ws_fast(x, y, mean_u, mean_v, pod_u, pod_v, dimensions, num_modes, vol_frac, bnd_idx, num_elem)
% Calculate spacing
hu = abs(mean(mean(x(1:end-1,:),2) - mean(x(2:end,:),2)));
hv = abs(mean(mean(y(:,1:end-1),1) - mean(y(:,2:end),1)));

% Add mean flows as mode zero
pod_u = [mean_u, pod_u];
pod_v = [mean_v, pod_v];

% TODO determine why these all need to be negative to match
pod_u = regroup(pod_u, dimensions);
pod_v = regroup(pod_v, dimensions);

pod_udx = zeros(size(pod_u));   % x derivative pod_u
pod_udy = zeros(size(pod_u));   % y derivative pod_u

pod_vdx = zeros(size(pod_v));   % x derivative pod_v
pod_vdy = zeros(size(pod_v));   % y derivative poy_v

for i = 1:size(pod_u,3)
   [pod_udy(:,:,i),pod_udx(:,:,i)] = gradient(pod_u(:,:,i), hu, hv);
   [pod_vdy(:,:,i),pod_vdx(:,:,i)] = gradient(pod_v(:,:,i), hu, hv);
end

pod_u = reshape(pod_u, [], num_modes);
pod_v = reshape(pod_v, [], num_modes);

% Convert all matrices to vectors for later Calculation
pod_udx     = reshape(pod_udx, num_elem, num_modes);
pod_udy     = reshape(pod_udy, num_elem, num_modes);
pod_vdx     = reshape(pod_vdx, num_elem, num_modes);
pod_vdy     = reshape(pod_vdy, num_elem, num_modes);

% Strip boundaries
[pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac] = ...
    strip_boundaries(bnd_idx, pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac);
end