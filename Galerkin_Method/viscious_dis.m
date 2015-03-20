function modal_eddy_vis = viscious_dis(modal_amp_flux, modal_amp_mean, num_modes, lambda2, l, q, Re0)
% Current Naive implementation of viscious disappation term, may want to
% explore other proposed viscous dissapation methods

modal_eddy_vis = zeros(num_modes, 1);
q_sum = zeros(num_modes,num_modes);
for i = 2:num_modes
    q_temp = reshape(q(i,:), num_modes, num_modes);
    for j = 1:num_modes;
        for k = 1:num_modes;
            q_sum(j,k) = q_temp(j, k)* mean(modal_amp_flux(:,j).*modal_amp_flux(:,k).*modal_amp_flux(:,i));
            if modal_amp_mean(j) ~= 0 && k == i
                q_sum(j,k) = q_sum(j,k) + q_temp(j, k)*mean(modal_amp_mean(j)*lambda2(i));
            elseif modal_amp_mean(k) ~= 0 && j == i
                q_sum(j,k) = q_sum(j,k) + q_temp(j, k)*mean(modal_amp_mean(k)*lambda2(i));
            end
        end
    end
    l_sum = l(i,i)*lambda2(i);
    modal_eddy_vis(i) = -(1/Re0 + sum(sum(q_sum))/sum(l_sum));
end
end