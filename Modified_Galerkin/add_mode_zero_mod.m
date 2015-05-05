function [modal_amp, pod_u, pod_v] = ...
    add_mode_zero_mod(modal_amp, pod_u, pod_v,  mean_u, mean_v)
% Add an additional zeroth mode corresponding to the flow mean

% Number of snapshots
num_images = size(modal_amp,1);

% Add mean and flucuating amplitudes
modal_amp = [ones(num_images,1), modal_amp];

% Add mode zero
pod_u = [mean_u, pod_u];
pod_v = [mean_v, pod_v];
end