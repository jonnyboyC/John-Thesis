function pod_vor = calc_pod_vor2(u, v, dimensions, cutoff)
% Calculate vorticity values using matlabs built in curl function
pod_vor = zeros(dimensions(1), dimensions(2), cutoff);
for i = 1:cutoff
    [pod_vor(:,:,i), ~] = curl(u(:,:,i), v(:,:,i));
end
end