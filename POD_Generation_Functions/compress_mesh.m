function [X, U] = compress_mesh(X, U)
% COMPRESS_MESH reduce a mesh by a factor of 2 in both directions
%
% [x, y, u, v] = COMPRESS_MESH(x, y, u, v) given raw PIV data provided by
% x, y, u, v average the quanitites in a 2x2 grid into a new mesh of halve
% the original resolution

[x, u] = flow_comps(X, U);
comps = flow_ncomps(X);
dims = flow_dims(X);

bnd_idx = flow_boundaries(U);
bnd_idx = double(bnd_idx <= 0);

ranges = zeros(1,dims);

% determine new range of image, i.e make even grid size 
for i = 1:comps
    ranges(i) = 2*floor(size(X.(x{i}),i)/2);
end

% Compressed dimensions
comp_dims = ranges./2;

% Preallocate masks fil with zeros
idx = cell(pow2(dims),1);
shifts = get_shifts(dims);

% Get shifted submatrices
for i = 1:length(idx)
    idx{i} = flow_index(repmat({[1, 2, -1]}, 1, dims), 1:dims, X.(x{1}));
    for j = 1:dims
        idx{i}{j} = idx{i}{j} + shifts(i,j);
    end
end

for i = 1:comps
    X_temp.(x{i}) = 0;
end
bnd_mask = false(comp_dims);

% Combine submatrices to produce compressed mesh
for i = 1:length(idx)
    bnd_mask = bnd_mask + logical(bnd_idx(idx{i}{:}));
    for j = 1:comps
        X_temp.(x{j}) = X_temp.(x{j}) + X.(x{j})(idx{i}{:});
    end
end

for i = 1:comps
    X.(x{i}) = X_temp.(x{i});
    U_temp.(u{i}) = 0;
end

dims2 = ndims(U.(u{1}));
images = size(U.(u{1}),dims2);
copy = ones(1, dims2);
copy(end) = images;
bnd_mask = logical(repmat(bnd_mask, copy));

for i = 1:length(idx)
    idx{i} = [idx{i}; {1:images}];
end

% Compress mesh for velocity data
for i = 1:length(idx)
    for j = 1:comps
        U_temp.(u{j}) = U_temp.(u{j}) + U.(u{j})(idx{i}{:});
    end
end

for i = 1:comps
   U_temp.(u{i})(bnd_mask) = 0;
   U.(u{i}) = U_temp.(u{i})/length(idx);
end

end