function da = system_odes2(~, a, model_coef, ni, modal_TKE, nonlinear)
% System of ode for time evolution of POD/Galerkin system
%
% SYSTEM_ODES(~, A, MODEL_COEF) returns the the changes in modal amplitudes
% for the system modes

num_modes = (size(model_coef,1));
da = zeros(num_modes, 1);

if nonlinear
    TKE = sum(1/2*a(2:end).^2);  
    ni = ni*sqrt(TKE/modal_TKE);
end

% Mean velocity mode remains at amplitude 1
da(1) = 0;
range = 2:num_modes;
idx = 1;

for j = 1:num_modes
    idx = idx + 1;
    da(range) = da(range) + model_coef(range, idx).*repmat(a(j),length(range),1).*ni(range);
end

% index upper half of matrix
for j = 1:num_modes
    for k = j:num_modes
        idx = idx + 1;
        da(range) = da(range) + model_coef(range,idx).*repmat(a(j)*a(k),length(range),1);
    end
end
