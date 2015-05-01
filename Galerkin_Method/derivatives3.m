function [dx, dy, d2x, d2y] = derivatives3(var, bnd_idx, bnd_x, bnd_y, x, y, dimensions)
% Find derivatives of variable var

% Figure out how many loops are needed
number2calc = size(var, 2);

% Prefill outputs
dx = zeros(dimensions(1), dimensions(2), number2calc);
dy = zeros(dimensions(1), dimensions(2), number2calc);
d2x = zeros(dimensions(1), dimensions(2), number2calc);
d2y = zeros(dimensions(1), dimensions(2), number2calc);

% ensure x and y and in 2D
x = reshape(x, dimensions);
y = reshape(y, dimensions);

% Select finite elemente method to be used at each grid point
[methodsX, methodsY] = select_method(bnd_idx, bnd_x, bnd_y, dimensions, true);

% Precalculate stencil for each point
[stencilX, stencilY] = generate_stencil(x, y, methodsX, methodsY, dimensions);

% Calculate derivative terms
for i = 1:number2calc
    % Padding to allow for simultaneous calculation of derivatives
    padded_var1 = [zeros(4,dimensions(2)); reshape(var(:,i), dimensions); zeros(4,dimensions(2))];
    padded_var2 = [zeros(dimensions(1),4), reshape(var(:,i), dimensions), zeros(dimensions(1),4)];
    
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