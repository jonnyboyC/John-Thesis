function [pod_u, pod_v, lambda2, eig_func] = calc_eig_modes2(co_var, num_modes, ...
    u_flux, v_flux)

co_var_gpu = gpuArray(co_var);
u_flux_gpu = gpuArray(u_flux);
v_flux_gpu = gpuArray(v_flux);
sz_co_var = length(co_var_gpu);
% Perform single value decomposition to get empirial eigenfunctions for
% first num_modes singular values
[eig_func, lambda2, qa] = svd(co_var_gpu);

% Normalize
eig_func = eig_func*realsqrt(lambda2)*sqrt(sz_co_var);

% Produce pod modes
pod_u = (u_flux_gpu*qa);
pod_v = (v_flux_gpu*qa);  

% Normalize pod modeses
pod_u = pod_u/realsqrt(sz_co_var*lambda2)';
pod_v = pod_v/realsqrt(sz_co_var*lambda2)';

pod_u = pod_u(:,1:num_modes);
pod_v = pod_v(:,1:num_modes);

% Derp 
eig_func = gather(eig_func);
lambda2 = gather(lambda2);
pod_u = gather(pod_u);
pod_v = gather(pod_v);

lambda2 = diag(lambda2);
end