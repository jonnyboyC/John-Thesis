function niu = viscious_dis(psi, num_modes, lambda2, l, q_dot, q)
% Current Naive implementation of viscious disappation term

% TODO what does nc mean, what is psi
nc = zeros(num_modes, num_modes);
niu = zeros(num_modes,1);
for k = 1:num_modes
    idx = 1;
    for i = 1:num_modes
        for j = 1:num_modes
            nc(i,j) = q(k,idx)*mean(psi(:,i).*psi(:,j).*psi(:,k));
            idx = idx+1;
        end
    end
    mf_e = q_dot(k,k)*lambda2(k);
    bc = l(k,k)*lambda2(k);
    qc = sum(sum(nc));
    uc = qc + mf_e;
    niu(k) = -uc/bc;
end