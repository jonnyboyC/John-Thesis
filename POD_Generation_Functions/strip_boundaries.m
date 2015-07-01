function varargout = strip_boundaries(bnd_idx, varargin)
% Check to ensure a cell is provided
if size(varargin,1) == 0
    error('Must provide at least one set of data');
end

% reshape boundary index to column vector
bnd_idx = reshape(bnd_idx, [], 1);

% determine number of points in active region
active = sum(bnd_idx == 1 | bnd_idx == 0);

% Strip boundaries and return results
varargout = cell(length(varargin), 1);
for i = 1:length(varargin)
    if ~ismatrix(varargin{i}) && ~isstruct(varargin{i})
        error('Provide Matrices must be 2D matrices');
    end
    
    % If struct assume one of the flow components
    if isstruct(varargin{i})
        varargout{i} = step_in(varargin{i}, bnd_idx, active);
    else
        num_images = size(varargin{i},2);
        mask = repmat(bnd_idx, 1, num_images);
        varargout{i} = varargin{i}(mask == 1 | mask == 0);
        varargout{i} = reshape(varargout{i}, active, num_images);
    end
end
end

% Allow repeated stepping into structures to strip inactive portions
function X = step_in(X, bnd_idx, active)
    % Get components
    comps = flow_ncomps(X);
    x = flow_comps(X);
    
    for i = 1:comps
        if isstruct(X.(x{i}))
            % If we get another structure step in further
            X.(x{i}) = step_in(X.(x{i}), bnd_idx, active);
        else
            % Else must be data and strip inactive points
            num_images = size(X.(x{i}),2);
            mask = repmat(bnd_idx, 1, num_images);
            X.(x{i}) = X.(x{i})(mask == 1 | mask == 0);
            X.(x{i}) = reshape(X.(x{i}), active, num_images);
        end
    end
end