function niu = viscious_dis_couplet(eig_func, num_modes, ...
                                    l_dot, l, q_2dot, q_dot, q, Re0)
num_images = size(eig_func, 1);
num_cutoff = size(q_2dot,1);    % Number of modes used for cutoff

% From Couplet et al. 2004
P_unres = zeros(num_modes, num_images);
D_res = zeros(num_modes, num_images);

% Calculate resolved terms for D
for k = 1:num_modes
    for i = 1:num_modes
        D_res(i,k)=l_dot(i);
        for j = 1:num_modes
            D_res(i,k)=D_res(i,k)+l(i,j)*eig_func(k,j);
        end
    end
end

% Calculate approximate unresolved terms for P using modes from num_modes
% to num_cutoff
for i = 1:num_modes
    q_temp=reshape(q(i,:),num_cutoff,num_cutoff);
    for j = num_modes+1:num_cutoff
        P_unres(i,:)=P_unres(i,:)+(q_dot(i,j)*eig_func(:,j)+l(i,j)*eig_func(:,j)/Re0)';
        for j1 = 1:j
            P_unres(i,:)=P_unres(i,:)+(q_temp(j,j1)*eig_func(:,j).*eig_func(:,j1))';
        end
    end
end

% niu = mean((P_unres'.*D_res),2)/mean((D_res^2),2);
niu2 = mean((eig_func(:,1:num_modes)'.*P_unres),2)./mean((eig_func(:,1:num_modes)'.*D_res),2);
% niu3 = mean((eig_func(:,1:num_modes)'.*P_unres*D_res),2)/mean((eig_func(:,1:num_modes)*D_res^2),2);

niu = niu2;
end