function [pod_u, pod_v, lambda2, eig_func, cutoff] = calc_eig_modes2(co_var, ...
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

lambda2 = diag(lambda2);

% Find the number of modes need to account for 99% of the energy
sum_mode_energy = cumsum(lambda2)./sum(lambda2);
cutoff = find((sum_mode_energy >= 0.99), 1);

% % Have to truncate for memory
% if cutoff > 600
%     cutoff = 600;
% end 

% Display cutoff energy content
fprintf('\nCutoff for Couplet Viscous Dissapation is %d mode at %3.4f % Energy\n\n', ...
    cutoff, sum_mode_energy(cutoff));

% Return truncated pod modes
pod_u = pod_u(:,1:cutoff);
pod_v = pod_v(:,1:cutoff);
end