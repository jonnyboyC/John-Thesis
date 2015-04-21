function [dx, d2x, dy, d2y] = derivatives3(var, bnd_idx, x, y, dimensions)
% Find derivatives of variable var

% Figure out how many loops are needed
number2calc = size(var, 2);

% Prefill
dx = zeros(dimensions(1), dimensions(2), number2calc);
dy = zeros(dimensions(1), dimensions(2), number2calc);
d2x = zeros(dimensions(1), dimensions(2), number2calc);
d2y = zeros(dimensions(1), dimensions(2), number2calc);

x = reshape(x, dimensions);
y = reshape(y, dimensions);

[methodsX, methodsY] = select_method(bnd_idx, dimensions);
[stencilX, stencilY] = generate_stencil(x, y, methodsX, methodsY, dimensions);

for i = 1:number2calc
    for j = 1:size(stencilX,3)
        dx(:,:,i) = dx(:,:,i) + reshape(var(:,i), dimensions).*stencilX(:,:,j);
        dy(:,:,i) = dy(:,:,i) + reshape(var(:,i), dimensions).*stencilY(:,:,j);
    end
    for j = 1:size(stencilX,3)
        d2x(:,:,i) = d2x(:,:,i) + dx(:,:,i).*stencilX(:,:,j);
        d2y(:,:,i) = d2y(:,:,i) + dy(:,:,i).*stencilY(:,:,j);
    end
end