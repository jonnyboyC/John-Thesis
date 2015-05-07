function da = system_odes_mod(~, a, model_coef)
% System of ode for time evolution of POD/Galerkin system
%
% SYSTEM_ODES(~, A, MODEL_COEF) returns the the changes in modal amplitudes
% for the system modes

% number of modes
num_modes = (size(model_coef,1));

% Contant terms
da = model_coef(:,1);

% Linear terms
da = da + sum(model_coef(:, 2:num_modes+1).*repmat(a',num_modes,1),2);

idx = num_modes+2;

% Quadractic terms, index upper half of matrix
for j = 1:num_modes
    for k = j:num_modes
        da = da + model_coef(:,idx).*repmat(a(j)*a(k),num_modes,1);
        idx = idx + 1;
    end
end