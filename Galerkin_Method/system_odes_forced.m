function da = system_odes_forced(t, a, model_coef, phase)
% SYSTEM_ODES_FORCED determine 1st derivative of the POD-Galerkin system fully
% vectorized with explict forcing
%
%   da = SYSTEM_ODES(t, A, MODEL_COEF, phase) 


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

% Forcing
da(2) = 0.06165*60/(2*pi)*cos(60/(2*pi)*t + phase);
da(3) = -0.06165*60/(2*pi)*sin(60/(2*pi)*t + phase);
end