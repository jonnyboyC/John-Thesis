function [bnd_X, bnd_idx] = refine_bounds(X, U, mean_U, direct, streamlines, ...
                                open_flow, update_bnds)
% REFINE_BOUNDS determine the location of various boundaries in the flow
%
% [bnd_x, bnd_y, bnd_idx] = REFINE_BOUNDS(x, y, u, v, mean_u, mean_v,
% direct, update_bnds) will launch a GUI interface to select the proper
% bounds if information is not already cached to a mat file

% Load boundaries if already determined
if exist([direct filesep 'Processed Data' filesep 'Mask' filesep 'Mask.mat'], 'file') && ~update_bnds
    load([direct filesep 'Processed Data' filesep 'Mask' filesep 'Mask.mat'], 'bnd_idx', 'bnd_X');
else
    % Determine flow boundaries
    bnd_idx = flow_boundaries(U, open_flow);
    
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
    [bnd_x, bnd_y, bnd_idx] = mask_gen(data, bnd_x, bnd_y, u, v, open_flow, streamlines);
    
    % Save values
    save([direct filesep 'Processed Data' filesep 'Mask' filesep 'Mask.mat'], 'bnd_idx', 'bnd_X', '-v7.3');
end

