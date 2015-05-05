function da = system_odes(~, a, model_coef)
% System of ode for time evolution of POD/Galerkin system
%
% SYSTEM_ODES(~, A, MODEL_COEF) returns the the changes in modal amplitudes
% for the system modes

num_modes = (size(model_coef,1));
da = zeros(num_modes, 1);

% Mean velocity mode remains at amplitude 1
range = 2:num_modes;

% Disspation terms
da(range) = sum(model_coef(range, 1:num_modes).*repmat(a(1:num_modes)',length(range),1),2);

% Convective terms, index upper half of matrix
idx = num_modes+1;

tic;
% index upper half of matrix
for j = 1:num_modes
    for k = j:num_modes
        da(range) = da(range) + model_coef(range,idx).*repmat(a(j)*a(k),length(range),1);
        idx = idx + 1;
    end
end
toc;

