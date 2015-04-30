function [C, L, Q, lambda, modal_amp] = term2order(l, q, vis, lambda, modal_amp)
% Convert from dissapative and convective form to interaction order
q = regroup_q(q);
l = repmat(vis, 1, size(l,1)).*l;
lambda = lambda(2:size(l,1));
modal_amp = modal_amp(:,2:end);

C = l(2:end,1) + q(2:end,1,1);
L = l(2:end, 2:end) + squeeze(q(1, 2:end, 2:end)) + squeeze(q(2:end, 1, 2:end));
Q = q(2:end, 2:end, 2:end);
end

function q = regroup_q(qi)
modes = size(qi,1);
q = zeros(modes, modes, modes);

for i = 1:modes
   q(:,:,i) = reshape(qi(i,:), [modes, modes]);
end
end