function bnd_idx = flow_boundaries(U)
% FLOW_BOUNDARIES determine flow boundaries
%
% bnd_idx = FLOW_BOUNDARIES(U) given raw flow data U, return
%   bnd_idx represented the different flow regions. 1 represents in flow, -1
%   represents uncaptured space, and 0 represented boundaries between the two

% Get information about components dimensions, and snapshots
u = flow_comps(U);
dims = flow_dims(U);
comps = flow_ncomps(U);
images = size(U.(u{1}), dims);

idx = flow_index({[1 1]}, dims(end), U);

% Intially set all to in flow
bnd_idx = zeros(size(squeeze(U.(u{1})(idx{:}))));

if ~ismatrix(bnd_idx)
    error(['Error currently bnd_idx has [only had a partial conversion to full 3d' ...
           'Currently can only handle 3d vector fields on a 2D plane']);
end

% TODO flow_boundaries are only partially converted still needs 2D plane as
% of now to work with edge detection

% Make all images with a pixel out of flow as part of boundary
for i = 1:images
    temp = zeros(size(bnd_idx));
    for j = 1:comps
        idx{end} = i;
        temp = temp + U.(u{j})(idx{:});
    end
    bnd_idx = bnd_idx + double(temp == 0);
end

% Set sets points with more one percent of image of zero velocity to out of
% the flow
bnd_idx(bnd_idx > images/100) = -1;
bnd_idx(bnd_idx >= 0) = 1;

% Uses Matlab community canny function for 2D-3D edge detection
bnd_edge = canny(bnd_idx, 0, 0, 'Centred', false);

% Use built in edge detection to boundary points change edge points to 0
bnd_idx(bnd_edge) = 0;

% Manual edge detection on image boundary
bnd_idx = manual_edge(bnd_idx);
end