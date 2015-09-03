function u_scale = u_scale_gen_shear(U, ~)
% U_SCALE_GEN_SHEAR calculate the free stream velocity of a shear layer
% assumed that flow is in the u direction 
%
%   u_scale = U_SCALE_GEN_SHEAR(U) 

dims = flow_dims(U);
u_scale = mean(U.('u'),dims);
u_scale = sort(abs(u_scale(:)));
u_scale = u_scale(floor(0.99*length(u_scale)));
end