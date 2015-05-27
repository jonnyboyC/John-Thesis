function[mean_u, mean_v, flux_u, flux_v] = clip_bounds(bnd_idx, mean_u, mean_v, flux_u, flux_v)
% CLIP_BOUNDS set all points inside the determined boundar to identically
% zero
%
% [mean_u, mean_v, flux_u, flux_v] = CLIP_BOUNDS(bnd_idx, mean_u, mean_v,
% flux_u, flux_v) set all points in mean in fluctuating components exactly
% to zero that have been determined to lie within the boundary

% bnd_idx == -1 represents values outside the flow
mean_u(bnd_idx <= 0) = 0;
mean_v(bnd_idx <= 0) = 0;

mask = repmat(bnd_idx, 1, 1, size(flux_u, 3));

flux_u(mask <= 0) = 0;
flux_v(mask <= 0) = 0;