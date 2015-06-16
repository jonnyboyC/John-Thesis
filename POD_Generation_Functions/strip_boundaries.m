function varargout = strip_boundaries(bnd_idx, num_images, varargin)
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
        comps = flow_ncomps(varargin{i});
        x = flow_comps(varargin{i});
        num_images = size(varargin{i}.(x{1}),2);
        mask = repmat(bnd_idx, 1, num_images);

        for j = 1:comps
            varargout{i}.(x{j}) = varargin{i}.(x{j})(mask == 1 | mask == 0);
            varargout{i}.(x{j}) = reshape(varargout{i}.(x{j}), active, num_images);
        end
    else
        num_images = size(varargin{i},2);
        mask = repmat(bnd_idx, 1, num_images);
        varargout{i} = varargin{i}(mask == 1 | mask == 0);
        varargout{i} = reshape(varargout{i}, active, num_images);
    end
end
end