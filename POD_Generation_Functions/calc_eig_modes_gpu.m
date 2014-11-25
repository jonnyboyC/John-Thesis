function [pod_u, pod_v, lambda2, eig_func] = calc_eig_modes_gpu(co_var, num_modes, ...
    u_flux, v_flux)

co_var_gpu = gpuArray(co_var);
sz_co_var = length(co_var_gpu);
% Perform single value decomposition to get empirial eigenfunctions for
% first num_modes singular values
[eig_func, lambda2, qa] = svd(co_var_gpu);
clear co_var_gpu

% Normalize
eig_func = eig_func*realsqrt(lambda2)*sqrt(sz_co_var);
norm_pod = realsqrt(sz_co_var*diag(lambda2));
norm_pod = norm_pod'*eye(size(lambda2,1))';

% Produce pod modes
u_flux_gpu = gpuArray(u_flux);
pod_u = (u_flux_gpu*qa);
clear u_flux_gpu
pod_u = pod_u/norm_pod;
pod_u = pod_u(:,1:num_modes);
pod_u = gather(pod_u);

% Normalize pod modeses
v_flux_gpu = gpuArray(v_flux);
pod_v = (v_flux_gpu*qa);  
clear v_flux_gpu
pod_v = pod_v/norm_pod;
pod_v = pod_v(:,1:num_modes);
pod_v = gather(pod_v);

% Derp 
eig_func = gather(eig_func);
lambda2 = gather(lambda2);
lambda2 = diag(lambda2);
end