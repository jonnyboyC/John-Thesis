function modal_eddy_vis = viscious_dis_couplet(modal_amp, num_modes, l, q, Re0)
num_images = size(modal_amp, 1);
num_cutoff = size(q,1);    % Number of modes used for cutoff

% From Couplet et al. 2004
P_unres = zeros(num_images, num_modes);
D_res = zeros(num_images, num_modes);

% Calculate resolved terms for D
sum_j = 1:num_modes;
for i = 2:num_modes
    D_res(:,i)= modal_amp(:,sum_j)*l(i,sum_j)';
end

% Calculate approximate unresolved terms for P using modes from num_modes
% to num_cutoff
k_sum = 1:num_cutoff;
for i = 2:num_modes
    q_temp = reshape(q(i,:), num_cutoff, num_cutoff);
    for j = num_modes+1:num_cutoff
        P_unres(:,i) = P_unres(:,i) + (l(i,j)*modal_amp(:,j)/Re0) ...
            + modal_amp(:,k_sum).*repmat(modal_amp(:,j),1,num_cutoff)*q_temp(j,k_sum)';
    end
end
sum_j = 2:num_modes;
modal_eddy_vis = zeros(num_modes, 1);

%modal_eddy_vis1 = mean((P_unres(:,sum_j).*D_res(:,sum_j)),1)'./mean((D_res(:,sum_j).^2),1)';
% modal_eddy_vis(2:end) = mean(((modal_amp(:,sum_j)).^2.*P_unres(:,sum_j).*D_res(:,sum_j)),1)'./mean((modal_amp(:,sum_j).*D_res(:,sum_j).^2),1)';
modal_eddy_vis(2:end) = mean((modal_amp(:,sum_j).*P_unres(:,sum_j)),1)'./mean((modal_amp(:,sum_j).*D_res(:,sum_j)),1)';
end