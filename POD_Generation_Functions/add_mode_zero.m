function [modal_amp, lambda, pod_u, pod_v, pod_vor] = ...
    add_mode_zero(modal_amp, lambda, pod_u, pod_v, pod_vor, mean_u, mean_v, mean_vor)
% ADD_MODE_ZERO add an additional mode corresponding to the mean flow, use
% edit for details on input, output

% Number of snapshots
num_images = size(modal_amp,2);

lambda = [0; lambda];

% Add mean and flucuating amplitudes
modal_amp = [ones(num_images,1), modal_amp];

% Add mode zero
pod_u = [mean_u, pod_u];
pod_v = [mean_v, pod_v];
pod_vor = [mean_vor, pod_vor];
end