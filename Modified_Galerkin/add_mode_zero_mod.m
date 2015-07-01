function [modal_amp, pod_U] = ...
    add_mode_zero_mod(modal_amp, pod_U,  mean_U)
% Add an additional zeroth mode corresponding to the flow mean

% Number of snapshots
num_images = size(modal_amp,1);

% Add mean and flucuating amplitudes
modal_amp = [ones(num_images,1), modal_amp];

comps = flow_comps(pod_U);
u = flow_comps(pod_U);

% Add mode zero
for i = 1:comps
    pod_U.(u{i}) = [mean_U.(u{i}), pod_U.(u{i})];
end

end