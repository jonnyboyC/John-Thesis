function dims = flow_dims(X)
% FLOW_DIMS the total number of dimensions in data
%
%   dims = FLOW_DIMS(X)

x = flow_comps(X);
dims = ndims(squeeze((X.(x{1}))));
end