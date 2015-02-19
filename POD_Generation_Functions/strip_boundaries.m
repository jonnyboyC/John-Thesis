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
varargout = cell(size(raw,2), 1);
for i = 1:size(raw,2)
    if length(size(raw{i})) == 3
        error('Provide Matrices must be 2D matrices');
    end
    num_images = size(raw{i},2);
    mask = repmat(bnd_idx, 1, num_images);
    varargout{i} = raw{i}(mask == 1 | mask == 0);
    varargout{i} = reshape(stripped{i}, active, num_images);
end
end