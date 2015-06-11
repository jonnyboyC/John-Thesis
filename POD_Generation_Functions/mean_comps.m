function varargout = mean_comps(varargin)
% MEAN_COMPS provides essentially same function as mean but for all
% components in X, U, etc
% 
%   MEAN_X = MEAN_COMPS(X) takes the mean along dimensions 1
%   MEAN_X = MEAN_COMPS(X, dim) takes mean along dimension dim
%   [MEAN_X, MEAN_U, ...] = MEAN_COMPS(X, U, ..., dim) take mean along 
%       several components on dimension dim

% check that proper inputs were provided
for i = 1:length(varargin)
    if ~isscalar(varargin{i}) || ~isstruct(varargin{i})
        error('Most provide component structs or dimension scalar');
    end
end

% Parse inputs into components and dim
if isscalar(varargin{end})
    dim = varargin{end};
    components = varargin{1:end-1};
else
    dim = 1;
    components = varargin{:};
end

% preallocate
varargout = cell(length(components,1));

% for each set of components provided
for i = 1:length(components)
    
    % Preallocate struct
    varargout{i} = ([]);
    
    % get field names and number of components
    x = flow_comps(components{i});
    comps = flow_ncomps(components{i});
    
    % take the mean
    for j = 1:comps
        varargout{i}.(x{j}) = mean(components.(x{j}), dim);
    end
end
