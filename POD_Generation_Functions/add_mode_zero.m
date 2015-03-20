function [modal_amp_mean, modal_amp_flux, lambda, pod_u, pod_v] = ...
    add_mode_zero(modal_amp_mean, modal_amp_flux, lambda, pod_u, pod_v,  mean_u, mean_v)
% Add an additional zeroth mode corresponding to the flow mean

% Number of snapshots
num_images = size(modal_amp_flux,2);

lambda = [0; lambda];

% Add mean and flucuating amplitudes
modal_amp_mean = [1, modal_amp_mean];
modal_amp_flux = [zeros(num_images,1), modal_amp_flux];

% Add mode zero
pod_u = [mean_u, pod_u];
pod_v = [mean_v, pod_v];
end