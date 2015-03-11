function modal_eddy_vis = viscious_dis_couplet(modal_amp, num_modes, l, q, Re0)
num_images = size(modal_amp, 1);
num_cutoff = size(q,1)+1;    % Number of modes used for cutoff
num_modes = num_modes + 1;

% From Couplet et al. 2004
P_unres = zeros(num_modes-1, num_images);
D_res = zeros(num_modes-1, num_images);

% Calculate resolved terms for D
for i = 1:num_modes-1
    for j = 1:num_modes
        D_res(i, :)=D_res(i, :) + l(i,j)*modal_amp(:,j)';
    end
end

% Calculate approximate unresolved terms for P using modes from num_modes
% to num_cutoff
for i = 1:num_modes-1
    q_temp = reshape(q(i,:), num_cutoff, num_cutoff);
    for j = num_modes+1:num_cutoff
        for k = 1:num_cutoff
        P_unres(i,:) = P_unres(i,:) + (l(i,j)*modal_amp(:,j)/Re0)' ...
       	             + (q_temp(j, k)*modal_amp(:,k).*modal_amp(:,j))';
        end
    end
end
sum_j = 2:num_modes;

% modal_eddy_vis1 = mean((P_unres.*D_res),2)./mean((D_res.^2),2);
% modal_eddy_vis3 = mean(((modal_amp(:,sum_j)').^2.*P_unres.*D_res),2)./mean((modal_amp(:,sum_j)'.*D_res.^2),2);

modal_eddy_vis = mean((modal_amp(:,sum_j)'.*P_unres),2)./mean((modal_amp(:,sum_j)'.*D_res),2);
end