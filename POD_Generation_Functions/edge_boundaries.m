function [bnd_x, bnd_y] = edge_boundaries(bnd_idx)
% EDGE_BOUNDARIES determine the velocity normal direction for the flow
%
% [bnd_x, bnd_y] = EDGE_BOUNDARIES(bnd_idx) given the boundary index found
% from FLOW_BOUNDARIES determine the velocity normal direction bnd_x, and
% bnd_y

% Determine the vector normal of the boundary
[bnd_y, bnd_x] = gradient(bnd_idx);

% Assume edges are also boundaries
bnd_y(:,1) = ones(size(bnd_idx,1),1);
bnd_y(:,end) = -ones(size(bnd_idx,1),1);

bnd_x(end,:) = -ones(size(bnd_idx,2),1);
bnd_x(1,:) = ones(size(bnd_idx,2),1);

% Set any values from gradient to either -1 or 1
bnd_y(bnd_y < 0) = -1;
bnd_y(bnd_y > 0) = 1;

bnd_x(bnd_x < 0) = -1;
bnd_x(bnd_x > 0) = 1;

% If a calculated value was in bnd_idx set to 0
bnd_y(bnd_idx == -1) = 0;
bnd_x(bnd_idx == -1) = 0;

exterior = false(size(bnd_idx));

bnd_x(bnd_idx == 1 & ~exterior) = 0;
exterior([1,size(bnd_idx,1)],:) = true;
exterior(:,[1,size(bnd_idx,2)]) = true;
bnd_y(bnd_idx == 1 & ~exterior) = 0;

end