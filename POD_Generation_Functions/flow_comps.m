function varargout = flow_comps(varargin)
% FLOW_COMPS get the field names of each flow variable
%
%   x = FLOW_COMPS(X)
%   u = FLOW_COMPS(U)
%   [x, u] = FLOW_COMPS(X, U)
%#ok<*AGROW>

varargout = cell(length(varargin));


for i = 1:length(varargin)
    temp = fieldnames(varargin{i});
    idx = strcmp(temp, 'direct');
    temp(idx) = [];
    varargout{i} = temp;
end
end