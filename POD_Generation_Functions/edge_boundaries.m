function [bnd_X] = edge_boundaries(bnd_idx, X)
% EDGE_BOUNDARIES determine the velocity normal direction for the flow
%
% [bnd_x, bnd_y] = EDGE_BOUNDARIES(bnd_idx) given the boundary index found
% from FLOW_BOUNDARIES determine the velocity normal direction bnd_x, and
% bnd_y

x = flow_comps_ns(X);

% Determine the vector normal of the boundary
switch ndims(bnd_idx) % TODO make this more abstract
    case 1
        bnd_X.(x{1}) = gradient(bnd_idx);
    case 2
        [bnd_X.(x{1}), bnd_X.(x{2})] = gradient(bnd_idx);
    case 3
        [bnd_X.(x{1}), bnd_X.(x{2}), bnd_X.(x{3})] = gradient(bnd_idx);
end

% Get the number of comps
comps = flow_ncomps(bnd_X);
dimensions = size(X.(x{1}));

% Assume edges are also boundaries
for i = 1:comps
    % Get for first and last index on a dimension
    
    idx_1 = flow_index([1, 1], i, bnd_X);
    idx_end = flow_index([0, 0], i, bnd_X);
    
    %TODO FIX FIX FIX
%     if i == 1
%         idx_1 = struct_index([1, 1], 2, bnd_X);
%         idx_end = struct_index([0, 0], 2, bnd_X);
%     else
%         idx_1 = struct_index([1, 1], 1, bnd_X);
%         idx_end = struct_index([0, 0], 1, bnd_X);
%     end
    
    % Get a dimensions - 1 set of ones and fill 
    dimensions_temp = dimensions;
    dimensions_temp(i) = 1;
    bnd_X.(x{i})(idx_1{:}) = ones(dimensions_temp);
    bnd_X.(x{i})(idx_end{:}) = -ones(dimensions_temp);
    
    % Set partial values to integars
    bnd_X.(x{i})(bnd_X.(x{i}) < 0) = -1;
    bnd_X.(x{i})(bnd_X.(x{i}) > 0) = 1;
    
    bnd_X.(x{i})(bnd_idx == -1) = 0;
end


% Generate a comps dimension matrix with true on the exterior
exterior = false(size(bnd_idx));

for i = 1:comps
    idx_1 = struct_index([1, 1], i, exterior);
    idx_end = struct_index([0, 0], i, exterior);
    
    % Set exterior of mesh to true
    exterior(idx_1{:}) = true;
    exterior(idx_end{:}) = true;
    
    % If radial grid revert any connected location
    if all(structfun(@(x) isequal(x(idx_1{:}), x(idx_end{:})), X))
        exterior(idx_1{:}) = false;
        exterior(idx_end{:}) = false;
    end
    if all(structfun(@(x) range(x(idx_1{:})) == 0, X))
        exterior(idx_1{:}) = false;
    end
    if all(structfun(@(x) range(x(idx_end{:})) == 0, X))
        exterior(idx_end{:}) = false;
    end
end

% Remove all points that are listed in the flow and not on the exterior
for i = 1:comps
    bnd_X.(x{i})(bnd_idx == 1 & ~exterior) = 0;
end

% Check to make sure we are not working with a polar coordinate field


end