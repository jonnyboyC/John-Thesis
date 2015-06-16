function[mean_U, flux_U] = clip_bounds(bnd_idx, mean_U, flux_U, copy)
% CLIP_BOUNDS set all points inside the determined boundar to identically
% zero
%
% [mean_u, mean_v, flux_u, flux_v] = CLIP_BOUNDS(bnd_idx, mean_u, mean_v,
% flux_u, flux_v) set all points in mean in fluctuating components exactly
% to zero that have been determined to lie within the boundary

% bnd_idx == -1 represents values outside the flow
comps = flow_ncomps(mean_U);
u = flow_comps(mean_U);
mask = repmat(bnd_idx, copy);

for i = 1:comps
    mean_U.(u{i})(bnd_idx <= 0) = 0;
    flux_U.(u{i})(mask <= 0) = 0;
end