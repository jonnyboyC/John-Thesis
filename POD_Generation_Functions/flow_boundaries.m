function bnd_idx = flow_boundaries(u, v, open_flow)
% FLOW_BOUNDARIES determine flow boundaries
%
% bnd_idx = FLOW_BOUNDARIES(u,v) given raw flow data u and v, return
% bnd_idx represented the different flow regions. 1 represents in flow, -1
% represents uncaptured space, and 0 represented boundaries between the two

% Intially set all to in flow
bnd_idx = ones(size(u(:,:,1)));

if open_flow
    return;
end

% Make all images with a pixel out of flow as part of boundary
for i = 1:size(u,3)
    bnd_idx = bnd_idx + double((u(:,:,i) + v(:,:,i) == 0));
end

% Set all points with at least one image not captured as in boundary
bnd_idx(bnd_idx ~= 1) = -1;

% Use built in edge detection to boundary points change edge points to 0
bnd_idx(edge(bnd_idx, 'canny')) = 0;

% Manual edge detection on image boundary
bnd_idx = manual_edge(bnd_idx);
end