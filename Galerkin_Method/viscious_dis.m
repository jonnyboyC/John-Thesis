function [modal_eddy_vis, global_eddy_vis] = viscious_dis(modal_amp, num_modes, lambda2, l, q, Re0)
% Current Naive implementation of viscious disappation term, may want to
% explore other proposed viscous dissapation methods

modal_eddy_vis = zeros(num_modes, 1);
q_globe = 0;
l_globe = 0;
q_sum = zeros(num_modes,num_modes);
for i = 2:num_modes
    q_temp = reshape(q(i,:), num_modes, num_modes);
    for j = 1:num_modes;
        for k = 1:num_modes;
            q_sum(j,k) = q_temp(j, k)* mean(modal_amp(:,j).*modal_amp(:,k).*modal_amp(:,i));
        end
    end
    l_sum = l(i,i)*lambda2(i);
    modal_eddy_vis(i) = -(1/Re0 + sum(sum(q_sum))/sum(l_sum));
    q_globe = q_globe + sum(sum(q_sum));
    l_globe = l_globe + sum(l_sum);
end
global_eddy_vis = repmat(-1/Re0 + q_globe/l_globe, num_modes,1);
global_eddy_vis(1) = 0;
end