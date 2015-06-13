function [X, U] = compress_mesh(X, U, open_flow)
% COMPRESS_MESH reduce a mesh by a factor of 2 in both directions
%
% [x, y, u, v] = COMPRESS_MESH(x, y, u, v) given raw PIV data provided by
% x, y, u, v average the quanitites in a 2x2 grid into a new mesh of halve
% the original resolution

[x, u] = flow_comps(x, u);
comps = flow_ncomps(X);
dims = flow_dims(X);

bnd_idx = flow_boundaries(U, open_flow);
bnd_idx = double(bnd_idx <= 0);

ranges = zeros(comps,1);

% determine new range of image, i.e make even grid size 
for i = 1:comps
    ranges(i) = 2*flor(size(X.(x{i}),i)/2);
end

% Compressed dimensions
comp_dims = ranges./2;

% Preallocate masks fil with zeros
masks = cell(pow2(dims), 2);
masks = {zeros(size(X.(x{1})))};

for i = 1:dims
    mask{i,1}(struct_index({[1, 2, -1],[2, 2, 0]},[1,2], X)) = 1;
    mask{i,1}(struct_index({[1, 2, -1],[2, 2, 0]},[1,2], X)) = 1;
end

mask1(1:2:range_x, 1:2:range_y) = 1;
mask2(2:2:range_x, 1:2:range_y) = 1;
mask3(1:2:range_x, 2:2:range_y) = 1;
mask4(2:2:range_x, 2:2:range_y) = 1;

mask1 = logical(mask1);
mask2 = logical(mask2);
mask3 = logical(mask3);
mask4 = logical(mask4);

% Create bounds mask, we want to element any new elements that had at least
% one element in the boundary previously
bnd_mask = logical(bnd_idx(mask1) + bnd_idx(mask2) + bnd_idx(mask3) + bnd_idx(mask4)); 
bnd_mask = reshape(bnd_mask, comp_dims);

% compress x and y
x = (x(mask1) + x(mask2) + x(mask3) + x(mask4))./4;
y = (y(mask1) + y(mask2) + y(mask3) + y(mask4))./4;

x = reshape(x, comp_dims);
y = reshape(y, comp_dims);

% Preallocate
u_temp = zeros([comp_dims, size(u,3)]);
v_temp = zeros([comp_dims, size(u,3)]);
for i = 1:size(u,3)
    ui = u(:,:,i);
    vi = v(:,:,i);
    
    % compress velocity data
    ui = (ui(mask1) + ui(mask2) + ui(mask3) + ui(mask4))./4;
    vi = (vi(mask1) + vi(mask2) + vi(mask3) + vi(mask4))./4;
    
    ui = reshape(ui, comp_dims);
    vi = reshape(vi, comp_dims);
    
    % remove boundary elements
    ui(bnd_mask) = 0;
    ui(bnd_mask) = 0;
    
    u_temp(:,:,i) = ui;
    v_temp(:,:,i) = vi;
end

u = u_temp;
v = v_temp;

end