function [pod_u, pod_v, lambda2, eig_func] = calc_eig_modes_gpu(co_var, num_modes, ...
    u_flux, v_flux)

gpuDevice(1);

co_var_gpu = gpuArray(co_var);
sz_co_var = length(co_var_gpu);
% Perform single value decomposition to get empirial eigenfunctions for
% first num_modes singular values
[eig_func, lambda2, qa] = svd(co_var_gpu);

% Normalize
eig_func = eig_func*realsqrt(lambda2)*sqrt(sz_co_var);
norm_pod = realsqrt(sz_co_var*lambda2);

norm_pod = gather(norm_pod);
eig_func = gather(eig_func);
lambda2  = gather(lambda2);
qa  = gather(qa);


% Produce pod modes
pod_u = (u_flux*qa);
pod_u = pod_u/norm_pod;
pod_u = pod_u(:,1:num_modes);


% Normalize pod modeses
pod_v = (v_flux*qa); 
pod_v = pod_v/norm_pod;
pod_v = pod_v(:,1:num_modes);

lambda2 = diag(lambda2);
end