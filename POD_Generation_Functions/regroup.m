function var_full = regroup(var, dimensions)
% REGROUP transform vector from of a variable back into full dimensions
%
%   var_full = REGROUP(var, dimensions)
images = size(var,2);

var_full = zeros([dimensions, images]);

var_idx = flow_index([1, 1], ndims(var_full), var_full);
% Convert from 2D to 3D
for i = 1:images
    var_idx{end} = i;
    var_full(var_idx{:}) = reshape(var(:,i), dimensions);
end

