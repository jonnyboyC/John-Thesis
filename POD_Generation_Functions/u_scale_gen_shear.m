function u_scale = u_scale_gen_shear(U, ~)
% Function to calculate u_scale
dims = flow_dims(U);
u_scale = mean(U.('u'),dims);
u_scale = sort(abs(u_scale(:)));
u_scale = u_scale(floor(0.99*length(u_scale)));
end