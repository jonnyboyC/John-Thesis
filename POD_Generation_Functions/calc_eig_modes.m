function [pod_u, pod_v, lambda2, eig_func] = calc_eig_modes(co_var, num_modes, ...
    u_data, v_data, mu_data, mv_data)

sz_co_var = length(co_var);

% Perform single value decomposition to get empirial eigenfunctions

[eig_func, lambda2, qa] = svds(co_var, num_modes);

% Normalize
eig_func = eig_func*sqrt(lambda2)*sqrt(sz_co_var);

for i = 1:size(u_data, 2);
    u_data(:,i) = u_data(:,i) - mu_data;
    v_data(:,i) = v_data(:,i) - mv_data;
end

pod_u = (u_data*qa(:,1:num_modes));
pod_v = (v_data*qa(:,1:num_modes));  

pod_u = pod_u(:,1:num_modes)/sqrt(sz_co_var*lambda2(1:num_modes,1:num_modes))';
pod_v = pod_v(:,1:num_modes)/sqrt(sz_co_var*lambda2(1:num_modes,1:num_modes))';

eig_func = eig_func(:,1:num_modes);
lambda2 = diag(lambda2);
end