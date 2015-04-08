function modal_eddy_vis = viscous_dis_couplet(modal_amp, num_modes, l, q, Re0)
% Determine modal visous disspation according the Couplet et al. 2004

num_images = size(modal_amp, 1);    % Number of snapshots
num_cutoff = size(q,1);             % Number of modes used for cutoff

% Prefill
Q_unres = zeros(num_images, num_modes); % Unresolved convective terms
L_unres = zeros(num_images, num_modes); % Unresolve dissapative terms
L_res   = zeros(num_images, num_modes); % Resolve dissapative terms

% Predefine sums
sum_i = 2:num_modes;
sum_j_res = 1:num_modes;
sum_j_unres = num_modes+1:num_cutoff;
k_sum = 1:num_cutoff;


% Calculate systems resolved dissapation
for i = sum_i
    L_res(:,i)   = L_res(:,i)   + modal_amp(:,sum_j_res)*l(i,sum_j_res)';
    L_unres(:,i) = L_unres(:,i) + modal_amp(:, sum_j_unres)*l(i,sum_j_unres)'/Re0;
end

% Determine approximately unresolved vicous and convective terms
for i = sum_i
    q_temp = reshape(q(i,:), num_cutoff, num_cutoff);
    for j = sum_j_unres
        Q_unres(:,i) = Q_unres(:,i) + modal_amp(:,k_sum).*repmat(modal_amp(:,j),1,num_cutoff)*q_temp(j,k_sum)';
    end
end

% Total unresolved system
P_unres = L_unres + Q_unres;

modal_eddy_vis = zeros(num_modes, 1);

% modal_eddy_vis1 = mean((P_unres(:,sum_i).*L_res(:,sum_i)),1)'./mean((L_res(:,sum_i).^2),1)';
% modal_eddy_vis(2:end) = mean(((modal_amp(:,sum_i)).^2.*P_unres(:,sum_i).*L_res(:,sum_i)),1)'./mean((modal_amp(:,sum_i).*L_res(:,sum_i)).^2,1)';
modal_eddy_vis(2:end) = (mean((modal_amp(:,sum_i).*P_unres(:,sum_i)),1)./mean((modal_amp(:,sum_i).*L_res(:,sum_i)),1))';
end