function da = system_odes2(~, a, model_coef)
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

% offset for convective terms 
idx = num_modes+1;

% Convective terms
da(range) = da(range) + ...
    squeeze(sum(sum(regroup(model_coef(range,idx:end)',[num_modes, num_modes])...
    .*repmat(a(1:num_modes)*a(1:num_modes)',1,1,length(range)),1),2));
end