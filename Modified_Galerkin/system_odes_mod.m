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

% offset for convective terms
idx = num_modes+2;

% Convective terms
da = da + ...
    squeeze(sum(sum(regroup(model_coef(:,idx:end)',[num_modes, num_modes])...
    .*repmat(a*a',1,1,num_modes),1),2));
end