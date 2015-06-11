function varargout = flow_ncomps(varargin)
% FLOW_NDIMS get the number of flow dimensions i.e x, y or x y z
%
%   xdims = FLOW_NDIMS(X)
%   udims = FLOW_NDIMS(U)
%   [xdims, udims] = FLOW_NDIMS(X,U)

varargout = cell(length(varargin));

for i = 1:length(varargin)
    varargout{i} = length(fieldnames(varargin{i})); 
end
end