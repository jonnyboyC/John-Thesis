function [s_X, nstencil] = generate_stencil_new(X, methods_X, dimensions)
% GENERATE_STENCIL generate the appropriate weights for each grid point
% when evaulating the derivtives
%
%   [s_X, nstencile] = GENERATE_STENCIL(X, methods_X, dimensions)
%       produces a matrix size = [dimesions, 9] representing the weights of
%       each point when included in derivatives calculations.


% Get structure components
xi = flow_comps_ip(X);
x_dims = flow_comps_ip(methods_X);

% Get flow dimensions
dims = flow_dims(X);

% TODO in the future could be dynamic
finite_order = 4;
nstencil = 2*finite_order + 1;

% dimension index
dim_idx = 1:dims;

% Generate stencile along each dimmension
for i = 1:dims
    
    % Prefill stencile structure
    temp_idx = [dimensions(circshift(dim_idx, i-1, 2)), 9];
    s_X.(x_dims{i}) = zeros(temp_idx);

    % Generate padding
    temp_dimensions = dimensions;
    temp_dimensions(i) = finite_order;
    padding = zeros(temp_dimensions);  
    
    % Pad mesh struture
    X.(xi{i}) = permute(cat(i, padding, X.(xi{i}), padding), circshift(dim_idx, -(i-1), 2)); 

    % Permute methods structure to have methods along dimension 1
    temp_methods = permute(methods_X.(x_dims{i}), circshift(dim_idx, i-1, 2));
    
    input = cell(1,nstencil);
    output = cell(1,nstencil);
    
    % Fill input cell
    for j = 1:length(input)
        idx = flow_index([j, j-9], 1, X.(xi{i}));
        input{j} = X.(xi{i})(idx{:});
    end
    
    [output{:}] = arrayfun(@stencil, input{:}, temp_methods);
    
    % Fill stencile structe from outputs
    idx = flow_index([1 1], ndims(s_X.(x_dims{i})), s_X.(x_dims{i}));
    for j = 1:length(output)
        idx{end} = j;
        s_X.(x_dims{i})(idx{:}) = output{j};
    end
    
    % Permute stencile back to original dims, remove NaN's
    s_X.(x_dims{i}) = -permute(s_X.(x_dims{i}), [circshift(dim_idx, i-1, 2),dims+1]);
    s_X.(x_dims{i})(isnan(s_X.(x_dims{i})) | isinf(s_X.(x_dims{i}))) = 0;
end
end

function [s1, s2, s3, s4, s5, s6, s7, s8, s9] = stencil(x1, x2, x3, x4, x5, x6, x7, x8, x9, method)
% STENCIL determine the coefficients for each points at the point of
% interest and 4 points shifted in both direction that much be summed to
% get the derivative of interest

switch method
    case 0
        % In boundary
        s1 = 0; s2 = 0; s3 = 0; s4 = 0; s5 = 0; s6 = 0; s7 = 0; s8 = 0; s9 = 0;
    case 1
        % Central 4th order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = fourthOrder(x3, x4, x5, x6, x7, 3);
    case 2 
        % Forward 1st order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = firstOrder(x5, x6, 1);
    case 3 
        % Forward 2nd order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = secondOrder(x5, x6, x7, 1);
    case 4
        % Forward 3rd order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = thirdOrder(x5, x6, x7, x8, 1);
    case 5
        % Forward 4th order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = fourthOrder(x5, x6, x7, x8, x9, 1);
    case 6
        % Forward biased 3rd order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = thirdOrder(x4, x5, x6, x7, 2);
    case 7
        % Forward biased 4th order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = fourthOrder(x4, x5, x6, x7, x8, 2);
    case 8
        % 1 point stencil
        s1 = 0; s2 = 0; s3 = 0; s4 = 0; s5 = 0; s6 = 0; s7 = 0; s8 = 0; s9 = 0;
    case 9
        % Central 2nd order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = secondOrder(x4, x5, x6, 2);
    case 10
        % Backwards biased 3rd order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = thirdOrder(x3, x4, x5, x6, 3);
    case 11
        % Backwards biased 4th order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = fourthOrder(x2, x3, x4, x5, x6, 4);
    case 12
        % Backwards 1st order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = firstOrder(x4, x5, 2);
    case 13
        % Backwards 2nd order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = secondOrder(x3, x4, x5, 3);
    case 14
        % Backwards 3rd order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = thirdOrder(x2, x3, x4, x5, 4);
    case 15
        % Backwards 4th order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = fourthOrder(x1, x2, x3, x4, x5, 5);
end
end