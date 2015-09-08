function prod = inner_prod(var1, var2, volume)
% INNER_PROD take the L2 inner produce between two sets of data
%
%   prod = inner_prod(var1, var2, volume), inner product of var1 and var2 
%       over a mesh with volume
prod = (var1'*(var2.*(repmat(volume, 1, size(var2,2)))))';
end