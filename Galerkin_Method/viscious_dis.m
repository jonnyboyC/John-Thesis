function modal_eddy_vis = viscious_dis(modal_amp_flux, modal_amp_mean, num_modes, l, q, Re0)
% Current Naive implementation of viscious disappation term, may want to
% explore other proposed viscous dissapation methods

modal_eddy_vis = zeros(num_modes, 1);
q_sum = zeros(num_modes,num_modes);
l_sum = zeros(num_modes,1);
for i = 1:num_modes
    q_temp = reshape(q(i,:), num_modes, num_modes);
    for j = 1:num_modes;
        l_sum(j) = l(i,j)*mean(modal_amp_flux(:,j).*modal_amp_flux(:,i));
        for k = 1:num_modes;
            q_sum(j,k) = q_temp(j, k)* mean(modal_amp_flux(:,j).*modal_amp_flux(:,k).*modal_amp_flux(:,i));
            if modal_amp_mean(j) ~= 0
                q_sum(j,k) = q_sum(j,k) + q_temp(j, k)*...
                    mean(modal_amp_mean(j).*modal_amp_flux(:,k).*modal_amp_flux(:,i));
            end
            if modal_amp_mean(k) ~= 0
                q_sum(j,k) = q_sum(j,k) + q_temp(j, k)*...
                    mean(modal_amp_mean(k).*modal_amp_flux(:,j).*modal_amp_flux(:,i));
            end
        end
    end
    modal_eddy_vis(i) = -(1/Re0 + sum(sum(q_sum))/sum(l_sum));
end
end