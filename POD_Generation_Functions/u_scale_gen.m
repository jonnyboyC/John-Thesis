function u_scale = u_scale_gen(u)
% Function to calculate u_scale
u_scale = mean(u,3);
u_scale = sort(abs(u_scale(:)));
u_scale = u_scale(floor(0.99*length(u_scale)));
end