function [modal_eddy_vis, global_eddy_vis] = eddy_vis_avg(modal_amp, num_modes, lambda2, l, q, vis)
% EDDY_VIS_AVG calculate eddy visocity in an averaged fashsion
%
%   modal_eddy_vis, global_eddy_vis = EDDY_VIS_AVG(modal_amp, num_modes,
%       lambda2, l, q, vis)

% Preallocate
modal_eddy_vis = zeros(num_modes, 1);
q_globe = 0;
l_globe = 0;
q_sum = zeros(num_modes,num_modes);

% 
for i = 2:num_modes
    q_temp = reshape(q(i,:), num_modes, num_modes);
    for j = 1:num_modes
        for k = 1:num_modes
            q_sum(j,k) = q_temp(j, k)* mean(modal_amp(:,j).*modal_amp(:,k).*modal_amp(:,i));
        end
    end
    l_sum = l(i,i)*lambda2(i);
    
    % Modal energy balance
    modal_eddy_vis(i) = -(vis + sum(sum(q_sum))/l_sum);
    
    % Sum convective terms and visocous terms
    q_globe = q_globe + sum(sum(q_sum));
    l_globe = l_globe + sum(l_sum);
end

% Calculate global terms
global_eddy_vis = repmat(-(vis + q_globe/l_globe), num_modes,1);
global_eddy_vis(1) = 0;
end