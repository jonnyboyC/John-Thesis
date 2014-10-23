load('D:\shear layer\PIVData\Data1\POD Data\POD.mat', 'pod_u1', 'dimensions', ...
    'bnd_idx');
z = ones(dimensions);
for i = 1:14
    [~, ~] = visder(reshape(pod_u1(:,1), dimensions), dimensions, z, bnd_idx);
    [~, ~] = visder2(reshape(pod_u1(:,1), dimensions), dimensions, z, bnd_idx);
end
clear all