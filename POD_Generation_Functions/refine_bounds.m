function [bnd_x, bnd_y] = refine_bounds(x, y, mean_u, bnd_idx, bnd_x, bnd_y, direct, new_mask)
cmax = max(max(abs(mean_u)));
data.x = x;
data.y = y;
data.bnd_idx = bnd_idx;
data.bnd_x = bnd_x;
data.bnd_y = bnd_y;
data.pod = mean_u/cmax;

if exist([direct filesep 'Processed Data' filesep 'Mask' filesep 'Mask.mat'], 'file') && ~new_mask
    load([direct filesep 'Processed Data' filesep 'Mask' filesep 'Mask.mat'], 'bnd_x', 'bnd_y');
else
    [bnd_x, bnd_y] = mask_gen(data, bnd_idx, bnd_x, bnd_y);
    save([direct filesep 'Processed Data' filesep 'Mask' filesep 'Mask.mat'], 'bnd_x', 'bnd_y', '-v7.3');
end
