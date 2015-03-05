function pod_vor = calc_pod_vor2(u, v, dimensions, cutoff)
%% TODO do this more efficiently 
pod_vor = zeros(dimensions(1), dimensions(2), cutoff);
for i = 1:cutoff
    [pod_vor(:,:,i), ~] = curl(u(:,:,i), v(:,:,i));
end
end