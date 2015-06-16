function covariance = cal_covariance_mat(flux_U, vol_frac, bnd_idx, num_images)
% CAL_COVARIANCE_MAT Calculate covariance matrix of ensemble data
%
% CAL_COVARIANCE_MAT(FLUX_U, FLUX_V, VOL_FRAC, BND_IDX) calculates
% covariance matrix of ensembles images based on the flucatating components
% of ensembles based on volume occupied by each elemented dicatated by
% VOL_FRAC. Reduces problem sized by stripping points from within the
% boundary from BND_IDX

% Remove data within boundaries
[flux_U, vol_frac] = strip_boundaries(bnd_idx, num_images, flux_U, vol_frac);
vol_frac = repmat(vol_frac, 1, num_images);

comps = flow_ncomps(flux_U);
u = flow_comps(flux_U);

% In Turbulence, Coherent Structures, Dynamical... this is 1/M*XX'*psi
covariance_partials = zeros(num_images);
for i = 1:comps
    covariance_partials = covariance_partials + (flux_U.(u{i}).*vol_frac)'*flux_U.(u{i});
end
covariance = covariance_partials/num_images;
end