function modal_eddy_vis = viscious_dis(modal_amp, num_modes, lambda2, l, q, Re0)
% Current Naive implementation of viscious disappation term, may want to
% explore other proposed viscous dissapation methods

modal_eddy_vis = zeros(num_modes, 1);
q_sum = zeros(num_modes,num_modes);
num_modes = num_modes + 1;
for i = 1:num_modes-1
    q_temp = reshape(q(i,:), num_modes, num_modes);
    for j = 1:num_modes;
        for k = 1:num_modes;
            q_sum(j,k) = q_temp(j, k)*mean(modal_amp(:,j).*modal_amp(:,k).*modal_amp(:,i));
        end
    end
    l_sum = l(i,i)*lambda2(i);
    modal_eddy_vis(i) = -(1/Re0 + sum(sum(q_sum))/l_sum);
end

end