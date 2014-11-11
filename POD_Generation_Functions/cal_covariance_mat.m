function covariance = cal_covariance_mat(flux_u, flux_v, vol_frac)
% May need to include n_zones back in

num_images = size(flux_u,2);
vol_frac = vol_frac*ones(1,num_images);

% In Turbulence, Coherent Structures, Dynamical... this is 1/M*XX'*psi
covariance = (((flux_u.*vol_frac)'*flux_u) + ((flux_v.*vol_frac)'*flux_u))/num_images;
end