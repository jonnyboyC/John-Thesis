function vorticity = calc_vorticity2(u, v, dimensions, cutoff)
%% TODO do this more efficiently 
vorticity = zeros(dimensions(1), dimensions(2), cutoff);
for i = 1:cutoff
    [vorticity(:,:,i), ~] = curl(u(:,:,i), v(:,:,i));
end
end