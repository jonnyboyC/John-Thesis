function bnd_idx = manual_edge(bnd_idx)
% MANUAL_EDGE determine flow boundaries not determined using built in edge
% detection
% bnd_idx = MANUAL_EDGE(bnd_idx) updates bnd_idx returned from MATLAB's
% edge function to update points along edge of image.

% Index the exterior ring
exterior = false(size(bnd_idx));
exterior([1,size(bnd_idx,1)],:) = true;
exterior(:,[1,size(bnd_idx,2)]) = true;

% Determine if the ring one from the edge is part of the flow
interior = false(size(bnd_idx));
interior(1,2:end)   = bnd_idx(2,2:end) == 1;
interior(2:end,1)   = bnd_idx(2:end,2) == 1;
interior(2:end,end) = bnd_idx(2:end,end-1) == 1;
interior(end,2:end) = bnd_idx(end-1,2:end) == 1;

% same for the corners
interior(1,1)       = bnd_idx(2,2) == 1; 
interior(1,end)     = bnd_idx(2,end-1) == 1;       
interior(end,1)     = bnd_idx(end-1,2) == 1;     
interior(end,end)   = bnd_idx(end-1, end-1) == 1; 

% canny won't detect an edge on image edge, manually do this
bnd_idx(exterior & interior & bnd_idx == -1) = 0;
end