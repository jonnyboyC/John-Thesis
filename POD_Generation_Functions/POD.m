function [pod_U, lambda, modal_amp, cutoff] = POD(co_var, flux_U)
% POD_2D perform Proper Orthogonal Decompostion on a data matrix decribed
% in Sorvich's POD method of snapshots
%
% [pod_u, pod_v, lambda, modal_amp, cutoff] = POD_2D(co_var, flux_u,
% flux_v) given the fluctuation components FLUX_U, FLUX_V, and covariance
% matrix CO_VAR return POD modes, POD_U, POD_V eigenvalues LAMBDA,
% modal_amplitudes MODAL_AMP, and modes to 99% CUTOFF

%% Calculate pod modes, pod lambda values, and the left eigenvector
num_images = size(co_var,1);

% Perform single value decomposition to get empirial eigenfunctions for
% first num_modes singular values
[modal_amp, lambda, ~] = svd(co_var);

% Calculate singular value of data matrix
sigma = sqrt(lambda*num_images);

comps = flow_ncomps(flux_U);
u = flow_comps(flux_U);

for i = 1:comps
    pod_U.(u{i}) = (flux_U.(u{i})*modal_amp) /sigma';
end

% Normalize
modal_amp = modal_amp*sigma;

% pull off diagonal
lambda = diag(lambda);

% Find the number of modes need to account for 99% of the energy
sum_mode_energy = cumsum(lambda)./sum(lambda);
cutoff = find((sum_mode_energy >= 0.99), 1);

% Display cutoff energy content
fprintf('%d modes need to capture %3.4f percent total energy\n\n', ...
    cutoff, sum_mode_energy(cutoff));

% Return truncated pod modes
for i = 1:comps
    pod_U.(u{i}) = pod_U.(u{i})(:,1:cutoff);
end
end