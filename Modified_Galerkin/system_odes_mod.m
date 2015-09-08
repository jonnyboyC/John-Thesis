function da = system_odes_mod(~, a, model_coef)
% SYSTEM_ODES determine 1st derivative of the POD-Galerkin system fully
% vectorized for the basis trasnformation method
%
%   da = SYSTEM_ODES_MOD(~, A, MODEL_COEF) 

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