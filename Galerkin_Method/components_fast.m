function [udx, udy, vdx, vdy, pod_udx, pod_udy, pod_vdx, pod_vdy, l_dot, l] ...
    = components_fast(x, y, mean_u, mean_v, pod_u, pod_v, dimensions, vol_frac, num_modes, num_elem)

hu = mean(mean(x(1:end-1,:),2) - mean(x(2:end,:),2));
hv = mean(mean(y(:,1:end-1),1) - mean(y(:,2:end),1));

% TODO determine why these all need to be negative to match
mean_u = reshape(-mean_u, dimensions(1), dimensions(2));
mean_v = reshape(-mean_v, dimensions(1), dimensions(2));
pod_u = regroup(-pod_u, dimensions);
pod_v = regroup(-pod_v, dimensions);

[udy,udx] = gradient(mean_u, hu, hv);   % gradient of mean u flow
[vdy,vdx] = gradient(mean_v, hu, hv);   % gradient of mean v flow
d2u = del2(mean_u, hu, hv)*4;           % laplacain of mean u flow
d2v = del2(mean_v, hu, hv)*4;           % laplacian of mean v flow

pod_udx = zeros(size(pod_u));   % x derivative pod_u
pod_udy = zeros(size(pod_u));   % y derivative pod_u
d2pod_u = zeros(size(pod_u));   % laplacain pod_u

pod_vdx = zeros(size(pod_v));   % x derivative pod_v
pod_vdy = zeros(size(pod_v));   % y derivative poy_v
d2pod_v = zeros(size(pod_v));   % laplacian pod_v

for i = 1:size(pod_u,3)
   [pod_udy(:,:,i),pod_udx(:,:,i)] = gradient(pod_u(:,:,i), hu, hv);
   [pod_vdy(:,:,i),pod_vdx(:,:,i)] = gradient(pod_v(:,:,i), hu, hv);
   d2pod_u(:,:,i) = del2(pod_u(:,:,i), hu, hv)*4;
   d2pod_v(:,:,i) = del2(pod_v(:,:,i), hu, hv)*4;
end

% Convert all matrix quantities to vectors for mean derivatives and
% laplacian
udx     = reshape(udx, num_elem, 1);
udy     = reshape(udy, num_elem, 1);
vdx     = reshape(vdx, num_elem, 1);
vdy     = reshape(vdy, num_elem, 1);
d2u     = reshape(d2u, num_elem, 1);
d2v     = reshape(d2v, num_elem, 1);

% Convert all matrix quantities to vectors for pod derivatives and
% laplacian
pod_udx     = reshape(pod_udx, num_elem, num_modes);
pod_udy     = reshape(pod_udy, num_elem, num_modes);
pod_vdx     = reshape(pod_vdx, num_elem, num_modes);
pod_vdy     = reshape(pod_vdy, num_elem, num_modes);
d2pod_u     = reshape(d2pod_u, num_elem, num_modes);
d2pod_v     = reshape(d2pod_v, num_elem, num_modes);

pod_u = reshape(pod_u, num_elem, num_modes);
pod_v = reshape(pod_v, num_elem, num_modes);

l_dot = inner_prod(d2u, pod_u, vol_frac) + inner_prod(d2v, pod_v, vol_frac);
l = inner_prod(d2pod_u, pod_u, vol_frac) + inner_prod(d2pod_v, pod_v, vol_frac);

end