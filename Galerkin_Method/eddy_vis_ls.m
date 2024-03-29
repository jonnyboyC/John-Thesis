function ls_eddy_vis = eddy_vis_ls(modal_amp, num_modes, modes, l, q, vis)
% EDDY_VIS_LS calculate eddy visocity in a least squares fashion according 
% to Couplet et a. 2004
%
%   ls_eddy_vis = EDDY_VIS_LS(modal_amp, num_modes, modes, l, q, vis)

num_images = size(modal_amp, 1);    % Number of snapshots
num_cutoff = size(q,1);             % Number of modes used for cutoff

% Prefill
Q_unres = zeros(num_images, num_modes); % Unresolved convective terms
L_unres = zeros(num_images, num_modes); % Unresolve dissapative terms
L_res   = zeros(num_images, num_modes); % Resolve dissapative terms

% Predefine sums
sum_i = modes;
sum_j_res = [1, modes];
sum_j_unres = 1:num_cutoff;
sum_j_unres(sum_j_res) = [];
k_sum = 1:num_cutoff;


% Calculate systems resolved dissapation
idx = 2;
for i = sum_i
    L_res(:,idx)   = L_res(:,idx)   + modal_amp(:,sum_j_res)*l(i,sum_j_res)';
    L_unres(:,idx) = L_unres(:,idx) + modal_amp(:, sum_j_unres)*vis*l(i,sum_j_unres)';
    idx = idx + 1;
end


% Determine approximately unresolved vicous and convective terms
idx = 2;
for i = sum_i
    q_temp = reshape(q(i,:), num_cutoff, num_cutoff);
    for j = sum_j_unres
        Q_unres(:,idx) = Q_unres(:,idx) + modal_amp(:,k_sum).*repmat(modal_amp(:,j),1,num_cutoff)*q_temp(j,k_sum)';
    end
    idx = idx + 1;
end

% Total unresolved system
P_unres = L_unres + Q_unres;

ls_eddy_vis = zeros(num_modes, 1);
ls_eddy_vis(2:end) = mean(((modal_amp(:,sum_i)).^2.*P_unres(:,2:num_modes).*L_res(:,2:num_modes)),1)'...
    ./mean((modal_amp(:,sum_i).*L_res(:,2:num_modes)).^2,1)';
end