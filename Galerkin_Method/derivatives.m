function [dX, d2X] = derivatives(U, bnd_idx, bnd_X, X, dimensions)
% DERIVATIVES take the first and 2nd derivatives of a set of variables
% using a nonuniform mesh 4th order finite difference method. DERIVATIVES
% selected the appropriate finite different method based on the boundaries
% defined in bnd_idx, bnd_x, and bnd_y
%
% [dx, dy] derivatives(var, bnd_idx, bnd_x, bnd_y, x, y, dimensions)
% calculate the first order derivatives
%
% [dx, dy, d2x, d2y] = derivative(var, bnd_idx, bnd_x, bnd_y, x, y,
% dimensions) calculate the first and 2nd order derivatives

% Figure out how many loops are needed
number2calc = size(U, 2);

comps = flow_ncomps(X);
[x, u] = flow_comp(X, U);

% Prefill outputs
for i = 1:comps
    dX.(x{i}) = zeros([dimensions, number2calc]);
    if nargout == 2
        d2X.(x{i}) = zeros([dimensions, number2calc]);
    end
end

for i = 1:comps
    X.(x{i}) = reshape(X.(x{i}), dimensions);
end

% Select finite elemente method to be used at each grid point
[methods_X] = select_method(bnd_idx, bnd_X, dimensions, true);

% Precalculate stencil for each point
[stencilX, stencilY] = generate_stencil(x, y, methodsX, methodsY, dimensions);

% Calculate derivative terms
for i = 1:number2calc
    % Padding to allow for simultaneous calculation of derivatives
    padded_var1 = [zeros(4,dimensions(2)); reshape(U(:,i), dimensions); zeros(4,dimensions(2))];
    padded_var2 = [zeros(dimensions(1),4), reshape(U(:,i), dimensions), zeros(dimensions(1),4)];
    
    % 1st order terms
    for j = 1:size(stencilX,3)
        dx(:,:,i) = dx(:,:,i) + padded_var1(10-j:end-j+1,:).*stencilX(:,:,j);
        dy(:,:,i) = dy(:,:,i) + padded_var2(:,10-j:end-j+1).*stencilY(:,:,j);
    end
    
    if nargout == 4
        % Padding to allow for simultaneous calculation of derivatives
        padded_dx = [zeros(4,dimensions(2)); dx(:,:,i); zeros(4,dimensions(2))];
        padded_dy = [zeros(dimensions(1),4), dy(:,:,i), zeros(dimensions(1),4)];

        % 2nd order terms
        for j = 1:size(stencilX,3)
            d2x(:,:,i) = d2x(:,:,i) + padded_dx(10-j:end-j+1,:).*stencilX(:,:,j);
            d2y(:,:,i) = d2y(:,:,i) + padded_dy(:,10-j:end-j+1).*stencilY(:,:,j);
        end
    end
end