function [bnd_x, bnd_y, bnd_idx] = refine_bounds(x, y, mean_u, mean_v, direct, update_bnds)
% Load boundaries if already determined
if exist([direct filesep 'Processed Data' filesep 'Mask' filesep 'Mask.mat'], 'file') && ~update_bnds
    load([direct filesep 'Processed Data' filesep 'Mask' filesep 'Mask.mat'], 'bnd_idx', 'bnd_x', 'bnd_y');
else
    % Determine flow boundaries
    bnd_idx = flow_boundaries(mean_u, mean_v);
    
    % Determine potential flow boundaries
    [bnd_x, bnd_y] = edge_boundaries(bnd_idx);
    
    % Set up data for plotting
    data.x = x;
    data.y = y;
    data.bnd_idx = bnd_idx;
    data.bnd_x = bnd_x;
    data.bnd_y = bnd_y;
    data.u = mean_u;
    data.v = mean_v;
    
    % Open GUI
    [bnd_x, bnd_y, bnd_idx] = mask_gen(data, bnd_x, bnd_y);
    
    % Save values
    save([direct filesep 'Processed Data' filesep 'Mask' filesep 'Mask.mat'], 'bnd_idx', 'bnd_x', 'bnd_y', '-v7.3');
end
