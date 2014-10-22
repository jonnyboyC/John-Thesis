function [dx, d2x, dy, d2y] = derivatives(var, dimensions, z, xxi, yxi, xet, yet, aj, bnd_idx)
% Figure out how many loops are needed
number2calc = size(var, 2);

% Prefill
dx = zeros(dimensions(1), dimensions(2), number2calc);
dy = zeros(dimensions(1), dimensions(2), number2calc);
d2x = zeros(dimensions(1), dimensions(2), number2calc);
d2y = zeros(dimensions(1), dimensions(2), number2calc);

for i = 1:number2calc
    
    [dxic, detc] = visder(reshape(var(:,i), size(z,1), size(z,2)),...
         dimensions(1), dimensions(2), z, bnd_idx);
     
    % Calculate approximation for first derivative
    [dx(:,:,i), dy(:,:,i)] = derivative_approx(dxic, detc, xxi, yxi, xet, yet, aj);

    % Calculate approximation for 2nd derivative for x
    [dxic, detc] = visder(dx(:,:,i), size(z,1), size(z,2), z, bnd_idx);
    [d2x(:,:,i), ~] = derivative_approx(dxic, detc, xxi, yxi, xet, yet, aj);

    % Calculate approximation for 2nd derivative for y
    [dxic, detc] = visder(dy(:,:,i), size(z,1), size(z,2), z, bnd_idx);
    [~, d2y(:,:,i)] = derivative_approx(dxic, detc, xxi, yxi, xet, yet, aj);
end
end

% Approximation of the derivative 
function [dx, dy] = derivative_approx(pxic, petc, xxi, yxi, xet, yet, aj)
    dx = aj.*(pxic.*yet - petc.*yxi);
    dy = aj.*(petc.*xxi - pxic.*xet);
end