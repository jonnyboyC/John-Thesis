function bnd_idx = flow_boundaries(U, open_flow)
% FLOW_BOUNDARIES determine flow boundaries
%
% bnd_idx = FLOW_BOUNDARIES(u,v) given raw flow data u and v, return
% bnd_idx represented the different flow regions. 1 represents in flow, -1
% represents uncaptured space, and 0 represented boundaries between the two

% Get information about components dimensions, and snapshots
u = flow_comps(U);
dims = flow_dims(U);
comps = flow_ncomps(U);
images = size(U.(u{1}), dims);

idx = struct_index({[1 1]}, dims(end), U);

% Intially set all to in flow
bnd_idx = ones(size(squeeze(U.(u{1})(idx{:}))));

if open_flow
    return;
end

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
bnd_idx(bnd_idx > 1) = 1;

% Use built in edge detection to boundary points change edge points to 0
bnd_idx(edge(bnd_idx, 'canny')) = 0;

% Manual edge detection on image boundary
bnd_idx = manual_edge(bnd_idx);
end