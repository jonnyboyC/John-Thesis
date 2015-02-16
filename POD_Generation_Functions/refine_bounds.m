function [bnd_x, bnd_y] = refine_bounds(x, y, mean_u, bnd_idx, bnd_x, bnd_y, direct, new_mask)
cmax = max(max(abs(mean_u)));
data.x = x;
data.y = y;
data.pod = mean_u/cmax;
data.cmax = 1;

if exist([direct '\Processed Data\Mask\Mask.mat'], 'file') && ~new_mask
    load([direct '\Processed Data\Mask\Mask.mat'], 'bnd_x', 'bnd_y');
else
    [bnd_x, bnd_y] = mask_gen(data, bnd_idx, bnd_x, bnd_y);
    save([direct '\Processed Data\Mask\Mask.mat'], 'bnd_x', 'bnd_y', '-v7.3');
end
