function niu = viscious_dis_couplet(modal_amp, num_modes, ...
                                    l_dot, l, q_2dot, q_dot, q, Re0, niu)
num_images = size(modal_amp, 1);
num_cutoff = size(q_2dot,1);    % Number of modes used for cutoff

% From Couplet et al. 2004
P_unres = zeros(num_modes, num_images);
D_res = zeros(num_modes, num_images);

% Calculate resolved terms for D
for k = 1:num_modes
    for i = 1:num_modes
        D_res(i,k)=l_dot(i);
        for j = 1:num_modes
            D_res(i,k)=D_res(i,k)+l(i,j)*modal_amp(k,j);
        end
    end
end

% Calculate approximate unresolved terms for P using modes from num_modes
% to num_cutoff
for i = 1:num_modes
    q_temp = reshape(q(i,:),num_cutoff,num_cutoff);
    for j = num_modes+1:num_cutoff
        P_unres(i,:) = P_unres(i,:)+(q_dot(i,j)*modal_amp(:,j)+l(i,j)*modal_amp(:,j)/Re0)';
        for k = 1:j
            P_unres(i,:) = P_unres(i,:)+(q_temp(j,k)*modal_amp(:,j).*modal_amp(:,k))';
        end
    end
end

niu = mean((P_unres.*D_res),2)./mean((D_res.^2),2);
niu2 = mean((modal_amp(:,1:num_modes)'.*P_unres),2)./mean((modal_amp(:,1:num_modes)'.*D_res),2);
niu3 = mean(((modal_amp(:,1:num_modes)').^2.*P_unres.*D_res),2)./mean((modal_amp(:,1:num_modes)'.*D_res.^2),2);

niu = niu2;
end