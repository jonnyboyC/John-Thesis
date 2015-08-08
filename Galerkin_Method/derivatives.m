function [UdX, Ud2X] = derivatives(U, bnd_idx, bnd_X, X, dimensions) 
% DERIVATIVES take the first and 2nd laplacian derivatives of a set of variables
% using a nonuniform mesh 4th order finite difference method. DERIVATIVES
% selectes the appropriate finite different method based on the boundaries
% defined in bnd_idx, bnd_X,
%
% UdX = derivatives(U, bnd_idx, bnd_X, X, X_direct, uniform, dimensions)
% calculate the first order derivatives
%
% [UdX, Ud2X] = derivative(U, bnd_idx, bnd_X, X, X_direct, uniform, dimensions)
% calculate the first and 2nd order derivatives
%
% TODO currently only calculates Lapalcian terms for 2nd order could extend
% to calculate a full Hessian fairly easily


%% Preprocessing 


% Get information on the components and the grid provided

[x, u] = flow_comps_ip(X, U);    % Get grid components i.e x, y etc.
dims = flow_dims(X);        % Number of active dimensions
full_dims = size(X.(x{1})); % Mesh dimensions


% TODO can only work with 2D grids at the moment
if dims > 2
    error(['This function can only handle two dimensional derivatives, a modification' ...
        'is needed to solve for 3D']);
end

% Determine number of vector fields present
fields = size(U.(u{1}), 2);

% Prefill velocity derivatives
for i = 1:dims
    UdX.(u{i}) = struct([]);
    for j = 1:dims
        x_temp.(x{j}) = zeros([dimensions, fields]); 
    end
    UdX.(u{i}) = x_temp;
    
    if nargout == 2
        Ud2X.(u{i}) = struct([]);
        for j = 1:dims
            x_temp.(x{j}) = zeros([dimensions, fields]);
        end
        Ud2X.(u{i}) = x_temp;
    end
end

% Convert to original dimensions if in a stacked vector
for i = 1:dims
    X.(x{i}) = reshape(X.(x{i}), dimensions);
end

%% Method selection and parametric transform setup

% Select finite element method and produce mesh stenciles
methods_X = select_method(bnd_idx, bnd_X, dimensions, x);

% Set parametric transform
    
% Copy dimensions for readability
xi = {'xi', 'eta', 'phi'};

% Create meshs of indices 
for i = 1:dims
    temp_dims = full_dims;

    setup = ones(1,dims);
    setup(i) = temp_dims(i);
    indices = zeros(setup);
    indices(:) = 1:temp_dims(i);
    temp_dims(i) = 1;

    Xi.(xi{i}) = repmat(indices, temp_dims);
end
xi = flow_comps(Xi);

% Select finite element method and produce mesh stenciles
[methods_Xi] = select_method(ones(dimensions), bnd_X, dimensions, xi);
[stencil_Xi, nstencil] = generate_stencil(Xi, methods_Xi, dimensions);

% Preallocate indices derivatives
for i = 1:dims
    XdXi.(xi{i}) = struct([]);
    for j = 1:dims
        x_temp.(x{j}) = zeros(dimensions);
    end
    XdXi.(xi{i}) = x_temp;
end

for i = 1:dims
    % Generate padding to vectorize derivative calculation
    temp_dimensions = dimensions;
    temp_dimensions(i) = floor(nstencil/2);
    padding = zeros(temp_dimensions);

    % index the stencil matrix
    stencil_idx = flow_index([1 1], ndims(stencil_Xi.(xi{1})), stencil_Xi);

    for j = 1:dims
        % Pad indices
        padded_var = cat(i, padding, reshape(X.(x{j}), dimensions), padding);

        % Calculte first order derivatives
        for k = 1:nstencil
            stencil_idx{end} = k;
            idx = flow_index([10-k, 1-k], i, padded_var);
            XdXi.(xi{j}).(x{i}) = XdXi.(xi{j}).(x{i}) + ...
                padded_var(idx{:}).*stencil_Xi.(xi{i})(stencil_idx{:});
        end
    end
end
% Get the final transform matrix taking inverse jacobian
transform = inverse_jacobian(XdXi, dimensions);

% Use computation domain mesh to generate stencil
[stencil_X, nstencil] = generate_stencil(Xi, methods_X, dimensions);

%% Calculate Derivatives

% Take the derivative in each dimension
for i = 1:dims
    
    % Generate padding to vectorize derivative calculation
    temp_dimensions = dimensions;
    temp_dimensions(i) = floor(nstencil/2);
    padding = zeros(temp_dimensions);
    
    % index the stencil matrix
    stencil_idx = flow_index([1 1], ndims(stencil_X.(x{1})), stencil_X);
    
    % Loop through each component and each field to be calculated
    for j = 1:dims
    for k = 1:fields

        % Pad variable
        padded_var = cat(i, padding, reshape(U.(u{j})(:,k), dimensions), padding);
        
        % 1st order terms
        for l = 1:nstencil
            stencil_idx{end} = l;
            idx = flow_index([10-l, 1-l], i, padded_var);
            UdX.(u{j}).(x{i})(:,:,k) = UdX.(u{j}).(x{i})(:,:,k) + ...
                    padded_var(idx{:}).*stencil_X.(x{i})(stencil_idx{:});
        end
        
        % Calculate 2nd order term if requested
        if nargout == 2
            
            % Padd variable derivative
            padded_var = cat(i, padding, UdX.(u{j}).(x{j})(:,:,k), padding);
            
            % 2nd order derivative
            for l = 1:nstencil
                stencil_idx{end} = l;
                idx = flow_index([10-l, 1-l], i, padded_var);
                Ud2X.(u{j}).(x{i})(:,:,k) = Ud2X.(u{j}).(x{i})(:,:,k) + ...
                        padded_var(idx{:}).*stencil_X.(x{i})(stencil_idx{:});
            end
        end
    end
    end
end

% Perform transform to get back to original domain
  
% Generate dimensiosn to copy matrices
copy = ones(1, dims+1);
copy(end) = fields;

% Copy matrix
UdX_copy = UdX;

% Combine terms to get transform
for i = 1:dims
    for j = 1:dims
        UdX.(u{i}).(x{j}) = 0;
        for k = 1:dims
            UdX.(u{i}).(x{j}) = UdX.(u{i}).(x{j}) + ...
                repmat(transform.(xi{k}).(x{j}), copy).*UdX_copy.(u{i}).(x{j});
        end
    end
end

if nargout == 2

    % Copy matrix
    Ud2X_copy = Ud2X;

    % Combine terms to get transform
    for i = 1:dims
        for j = 1:dims
            Ud2X.(u{i}).(x{j}) = 0;
            for k = 1:dims
                Ud2X.(u{i}).(x{j}) = Ud2X.(u{i}).(x{j}) + ...
                    repmat(transform.(xi{k}).(x{j}), copy).*Ud2X_copy.(u{i}).(x{j});
            end
        end
    end
end
end
