function covariance = cal_covariance_mat(u_data, v_data, mu_data, mv_data, vol_frac)
% May need to include n_zones back in

num_images = size(u_data,2);

%Remove mean velocity from snapshots 
u_data = u_data-mu_data*ones(1,num_images);
v_data = v_data-mv_data*ones(1,num_images);
vol_frac = vol_frac*ones(1,num_images);

% In Turbulence, Coherent Structures, Dynamical... this is 1/M*XX'*psi
covariance = (((u_data.*vol_frac)'*u_data) + ((v_data.*vol_frac)'*v_data))/num_images;
end