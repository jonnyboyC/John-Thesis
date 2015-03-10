function covariance = cal_covariance_mat2(flux_u, flux_v, vol_frac, bnd_idx)
%% Calculate covariance matrix of ensemble data, 
num_images = size(flux_u,2);

% Remove data within boundaries
[flux_u, flux_v, vol_frac] = strip_boundaries(bnd_idx, flux_u, flux_v, vol_frac);
vol_frac = vol_frac*ones(1,num_images);

% In Turbulence, Coherent Structures, Dynamical... this is 1/M*XX'*psi
covariance = (((flux_u.*vol_frac)'*flux_u) + ((flux_v.*vol_frac)'*flux_v))/num_images;
end