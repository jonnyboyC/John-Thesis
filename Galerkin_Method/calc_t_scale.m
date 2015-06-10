function [t_scale, tspan] = calc_t_scale(u_scale, l_scale, non_dim, tspan)
% CALC_T_SCALE determine the appropraite time scaling based on the velocity
% and length scales if the data was non-dimensionalized
%
%   [t_scale, tspan] = CALC_T_SCALE(u_scale, l_scale, non_dim, tspan);

if non_dim
    t_scale = u_scale/l_scale;
else
    t_scale = 1;
end
tspan = tspan*t_scale;

end