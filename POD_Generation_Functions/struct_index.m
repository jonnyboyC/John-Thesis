function varargout = struct_index(indices, dims, varargin)
% STRUCT_INDEX Index the structured values X, U, and POD
%
%   indices = STRUCT_INDEX(indices, dims, X) 
%   [indices1, indices2] = STRUCT_INDEX(indices, dims, X, POD)
%
%   Provide a cell array of indices such as {[1 10], [10 20]} and
%   a matrix of dims such as [1, 2]. Here the this would produce an index
%   equailvent to X(1:10, 10:20, :) where all other dimensions are taken to
%   be the : operator
%
%   To index similar to end use negative numbers where 0 = end and -2 =
%   end - 1 and so on
%
%   Return indices are returned as a cell array such that
%   X.(x(1))(indices{:}) provides the index

% Check that a cell was given for indices
if ~iscell(indices)
    error('Provide cell array for struct_index to index see help for details');
end

% Check that a vector was given for dims
if ~isvector(dims)
    error('Provide a row vector for dims see help for details');
end

% Check that a struct cooresponding to X, U, or POD
if any(~cellfun(@isstruct, varargin))
    error('Supply one of the component structure X, U, or POD');
end

% get data in correct order
[dims, idx] = sort(dims);
indices = indices(idx);

varargout = cell(length(varargin),1);

% Go through each variable provided
for i = 1:length(varargin)
    
    % get structs field names and dimensions
    x = flow_comps(varargin{i});
    full_dims = size(varargin{i}.(x{1}));
    varargout{i} = cell(length(full_dims), 1); 
    
    % Loop through each dimension where neglected dimensions are given :
    for j = 1:length(full_dims)
        
        % Check if index for dimensions was requested
        idx = find(dims == j);
        if ~isempty(idx)
            
            % Create equivalent to end statements or end - n statements
            if indices{idx}(1) <= 0
                indices{idx}(1) = full_dims(j) + indices{idx}(1);
            end
            if indices{idx}(2) <= 0
                indices{idx}(2) = full_dims(j) + indices{idx}(2);
            end
            
            % Check that we're actually in the dimension bounds
            if indices{idx}(1) < 0 || indices{idx}(1) > full_dims(j)
                error(['You have provided a set of indices out of bounds in indices %d' ...
                    'where you provided %d value with dimension length %d'], idx, indices{idx}(1), full_dims(j));
            end
            if indices{idx}(2) < 0 || indices{idx}(2) > full_dims(j)
                error(['You have provided a set of indices out of bounds in indices %d' ...
                    'where you provided %d value with dimension length %d'], idx, indices{idx}(1), full_dims(j));
            end
            
            % Add dimensions
            varargout{i}{j} = indices{idx}(1):indices{idx}(2);
        else
            % Equivalent to 1:end
            varargout{i}{j} = 1:full_dims(j);
        end
    end
end
end