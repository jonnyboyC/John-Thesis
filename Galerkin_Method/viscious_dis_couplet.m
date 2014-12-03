function niu = viscious_dis_couplet(eig_func, num_modes, lambda2,...
                                    l_dot, l, q_2dot, q_dot, q)
num_truncated = size(eig_func, 1);
num_unresolved = szie(q_dot,1);

for i = 1:num_unresolved
    
end

% From Couplet et al. 2004
P = zeroes(num_modes, num_truncated);
d = zeros(num_modes, num_truncated);

% Calculate resolved terms for d
for i = 1:num_truncated
    for j = 1:num_modes
        d(i,j) = l_dot(i) + l(i,1:num_modes).*eig_func(i,1:num_modes);
    end
end

% Calculate approximate unresolved terms for p using more nodes
for i = 1:num_truncated
    for j = 1:num_modes
        q_temp = reshape(q(j,:), num_unresolved, num_unresolved);
        for k = num_modes+1:num_unresolved
            
        end
    end
end



end