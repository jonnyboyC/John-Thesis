function [dx, d2x, dy, d2y] = derivatives2(var, x, y, dimensions)
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

for i = 1:number2calc     
    % Calculate approximation for first derivative
    dx(:,:,i) = DGradient(reshape(var(:,i), dimensions), x(:,1), 1, '2ndOrder');
    dy(:,:,i) = DGradient(reshape(var(:,i), dimensions), y(1,:), 2, '2ndOrder');

    % Calculate approximation for second derivative
    d2x(:,:,i) = DGradient(dx(:,:,i), x(:,1), 1, '2ndOrder');
    d2y(:,:,i) = DGradient(dy(:,:,i), y(1,:), 2, '2ndOrder');
end
end