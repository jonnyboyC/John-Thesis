function [dx, d2x, dy, d2y] = derivatives(var, dimensions, z, xxi, yxi, xet, yet, aj, bnd_idx)
% Find derivatives of variable var

% Figure out how many loops are needed
number2calc = size(var, 2);

% Prefill
dx = zeros(dimensions(1), dimensions(2), number2calc);
dy = zeros(dimensions(1), dimensions(2), number2calc);
d2x = zeros(dimensions(1), dimensions(2), number2calc);
d2y = zeros(dimensions(1), dimensions(2), number2calc);

parfor i = 1:number2calc
    
    % TODO look into parrallelizing visder
    [dxic, detc] = visder2(reshape(var(:,i), size(z,1), size(z,2)),...
         dimensions(1), dimensions(2), z, bnd_idx);
     
    % Calculate approximation for first derivative
    [dx(:,:,i), dy(:,:,i)] = derivative_approx(dxic, detc, xxi, yxi, xet, yet, aj);

    % Calculate approximation for 2nd derivative for x
    [dxic, detc] = visder2(dx(:,:,i), size(z,1), size(z,2), z, bnd_idx);
    [d2x(:,:,i), ~] = derivative_approx(dxic, detc, xxi, yxi, xet, yet, aj);

    % Calculate approximation for 2nd derivative for y
    [dxic, detc] = visder2(dy(:,:,i), size(z,1), size(z,2), z, bnd_idx);
    [~, d2y(:,:,i)] = derivative_approx(dxic, detc, xxi, yxi, xet, yet, aj);
end
end

% Approximation of the derivative 
function [dx, dy] = derivative_approx(dxic, detc, xxi, yxi, xet, yet, aj)
    dx = aj.*(dxic.*yet - detc.*yxi);
    dy = aj.*(detc.*xxi - dxic.*xet);
end