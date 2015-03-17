function [pod_u, pod_v, lambda2, modal_amp_mean, modal_amp_flux, cutoff] = ...
    calc_eig_modes2(co_var, flux_u, flux_v, mean_u, mean_v)
%% Calculate pod modes, pod lambda values, and the left eigenvector
num_images = size(co_var,1);

% Perform single value decomposition to get empirial eigenfunctions for
% first num_modes singular values
[modal_amp, lambda2, ~] = svd(co_var);

% Calculate singular value of data matrix
sigma = sqrt(lambda2*num_images);

% Produce pod modes
pod_u = (flux_u*modal_amp)/sigma';
pod_v = (flux_v*modal_amp)/sigma';  

% Normalize
modal_amp = modal_amp*sigma;

% Remove numerical errors
modal_amp_mean = mean(modal_amp, 1);
modal_amp_mean(modal_amp_mean < 1e-15) = 0;

% Decompose into mean and flucuating amplitudes
modal_amp_flux = zeros(size(modal_amp));
for i = 1:num_images
   modal_amp_flux(:,i) = modal_amp(:,i) - modal_amp_mean(i); 
end

% pull off diagonal
lambda2 = diag(lambda2);

% Find the number of modes need to account for 99% of the energy
sum_mode_energy = cumsum(lambda2)./sum(lambda2);
cutoff = find((sum_mode_energy >= 0.99), 1);

% Display cutoff energy content
fprintf('\nCutoff for Couplet Viscous Dissapation is %d mode at %3.4f percent total energy\n\n', ...
    cutoff, sum_mode_energy(cutoff));

% Return truncated pod modes
pod_u = pod_u(:,1:cutoff);
pod_v = pod_v(:,1:cutoff);
end