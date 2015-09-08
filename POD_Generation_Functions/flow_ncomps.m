function varargout = flow_ncomps(varargin)
% FLOW_NCOMPS get the number of flow dimensions i.e x, y or x y z
%
%   xdims = FLOW_NCOMPS(X)
%   udims = FLOW_NCOMPS(U)
%   [xdims, udims] = FLOW_NCOMPS(X,U)

varargout = cell(length(varargin));

for i = 1:length(varargin)
    if isfield(varargin{i}, 'direct');
        varargout{i} = length(fieldnames(rmfield(varargin{i}, 'direct'))); 
    else
        varargout{i} = length(fieldnames(varargin{i}));
    end
end
end