function da = system_odes_forced_NL(t, a, model_coef, niu, vis, modal_TKE, phase)
% SYSTEM_ODES_FORCED_NL determine 1st derivative of the POD-Galerkin system fully
% vectorized with explict forcing and non-linear eddy visocity
%
%   da = SYSTEM_ODES_NL(~, A, MODEL_COEF)

num_modes = (size(model_coef,1));
da = zeros(num_modes, 1);

TKE = sum(1/2*a(2:end).^2);  
niu = niu*sqrt(TKE/modal_TKE);
vis_total = niu + vis;

% Mean velocity mode remains at amplitude 1
range = 2:num_modes;

da(range) = sum(model_coef(range, 1:num_modes).*repmat(a(1:num_modes)',length(range),1),2)...
    .*vis_total(range);

idx = num_modes+1;

% index upper half of matrix
da(range) = da(range) + ...
    squeeze(sum(sum(regroup(model_coef(range,idx:end)',[num_modes, num_modes])...
    .*repmat(a(1:num_modes)*a(1:num_modes)',1,1,length(range)),1),2));

da(2) = 0.06165*60/(2*pi)*cos(60/(2*pi)*t + phase);
da(3) = -0.06165*60/(2*pi)*sin(60/(2*pi)*t + phase);
end
