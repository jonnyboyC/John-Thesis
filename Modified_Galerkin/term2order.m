function [C, L, Q] = term2order(l, q)
% Convert from dissapative and convective form to interaction order
C = l(2:end,1) + q(2:end,1,1);
L = l(2:end, 2:end) + q(1, 2:end, 2:end) + q(2:end, 1, 2:end);
Q = q(2:end, 2:end, 2:end);
end