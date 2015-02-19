function covariance = cal_covariance_mat2(flux_u, flux_v, vol_frac, bnd_idx)
%% Calculate covariance matrix of ensemble data, 
num_images = size(flux_u,2);

% Remove data within boundaries
stripped = strip_boundaries({flux_u, flux_v, vol_frac}, bnd_idx);

flux_u = stripped{1};
flux_v = stripped{2};
vol_frac = stripped{3};

% In Turbulence, Coherent Structures, Dynamical... this is 1/M*XX'*psi
covariance = (((flux_u.*vol_frac)'*flux_u) + ((flux_v.*vol_frac)'*flux_v))/num_images;

end