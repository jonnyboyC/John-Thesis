function vis = calc_vis(u_scale, l_scale, Re0, non_dim)
% CALC_VIS determine the flows visocity based on the Reynolds numbers,
% characterist velocity and length scale and if the results were
% non_dimensionalized
%
%   vis = CALC_VIS(u_scale, l_scale, Re0, non_dim)

if non_dim
    vis = 1/Re0;
else
    vis = (u_scale*l_scale)/Re0;
end 
end