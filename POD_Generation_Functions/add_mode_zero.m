function [modal_amp, lambda, pod_U] = ... , pod_VOR] = ...
    add_mode_zero(modal_amp, lambda, pod_U, mean_U) %pod_VOR, mean_VOR)
% ADD_MODE_ZERO add an additional mode corresponding to the mean flow, use
% edit for details on input, output

% Number of snapshots
num_images = size(modal_amp,2);

lambda = [0; lambda];

% Add mean and flucuating amplitudes
modal_amp = [ones(num_images,1), modal_amp];

comps = flow_ncomps(pod_U);
u = flow_comps(pod_U);

% Add mode zero
for i = 1:comps
    pod_U.(u{i}) = [mean_U.(u{i}), pod_U.(u{i})];
end
%pod_vor = [mean_vor, pod_vor];
end