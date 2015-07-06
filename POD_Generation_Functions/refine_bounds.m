function [bnd_X, bnd_idx] = refine_bounds(X, U, mean_U, direct, streamlines, update_bnds)
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
    bnd_idx = flow_boundaries(U);
    
    % Determine if boundaries are visual or real
    bnd_X = edge_boundaries(bnd_idx, X);
    
    % Set up data for plotting
    data.X = X;
    data.bnd_idx = bnd_idx;
    data.bnd_X = bnd_X;
    data.U = mean_U;
    
    % Open GUI
    [bnd_X, bnd_idx] = mask_gen(data, U, streamlines);
    
    % Save values
    save([direct filesep 'Processed Data' filesep 'Mask' filesep 'Mask.mat'], 'bnd_idx', 'bnd_X', '-v7.3');
end

