function [pod_u, pod_v, lambda2, eig_func] = calc_eig_modes2(co_var, num_modes, ...
    u_flux, v_flux)
%% Calculate pod modes, pod lambda values, and the left eigenvector

sz_co_var = length(co_var);

% Perform single value decomposition to get empirial eigenfunctions for
% first num_modes singular values
[eig_func, lambda2, qa] = svd(co_var);

% Normalize
eig_func = eig_func*realsqrt(lambda2)*sqrt(sz_co_var);

% Produce pod modes
pod_u = (u_flux*qa);
pod_v = (v_flux*qa);  

% Normalize pod modeses
pod_u = pod_u/sqrt(sz_co_var*lambda2)';
pod_v = pod_v/sqrt(sz_co_var*lambda2)';

% Return truncated pod modes
pod_u = pod_u(:,1:num_modes);
pod_v = pod_v(:,1:num_modes);

lambda2 = diag(lambda2);
end