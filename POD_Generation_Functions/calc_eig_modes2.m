function [pod_u, pod_v, lambda2, eig_func] = calc_eig_modes2(co_var, num_modes, ...
    u_data, v_data, mu_data, mv_data)

sz_co_var = length(co_var);

% Perform single value decomposition to get empirial eigenfunctions for
% first num_modes singular values
[eig_func, lambda2, qa] = svds(co_var, num_modes);

% Normalize
eig_func = eig_func*sqrt(lambda2)*sqrt(sz_co_var);

% Get velocity fluxations from instaneous velocity
for i = 1:size(u_data, 2);
    u_data(:,i) = u_data(:,i) - mu_data;
    v_data(:,i) = v_data(:,i) - mv_data;
end

% Produce pod modes
pod_u = (u_data*qa);
pod_v = (v_data*qa);  

% Normalize pod modeses
pod_u = pod_u/sqrt(sz_co_var*lambda2)';
pod_v = pod_v/sqrt(sz_co_var*lambda2)';

lambda2 = diag(lambda2);
end