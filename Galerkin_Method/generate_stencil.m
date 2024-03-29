function [s_X, nstencil] = generate_stencil(X, methods_X, dimensions)
% GENERATE_STENCIL generate the appropriate weights for each grid point
% when evaulating the derivtives
%
%   [s_X, nstencile] = GENERATE_STENCIL(X, methods_X, dimensions)
%       produces a matrix size = [dimesions, 9] representing the weights of
%       each point when included in derivatives calculations.

xi = flow_comps_ip(X);
x_dims = flow_comps_ip(methods_X);
x = X.(xi{1});
y = X.(xi{2});

methodsX = methods_X.(x_dims{1});
methodsY = methods_X.(x_dims{2});

sX = zeros([dimensions, 9]);
sY = zeros([dimensions(2), dimensions(1), 9]);

% Add padding
x = [zeros(4,dimensions(2)); x; zeros(4,dimensions(2))];
y = [zeros(dimensions(1),4), y, zeros(dimensions(1),4)]';
methodsY = methodsY';

[sX(:,:,1), sX(:,:,2), sX(:,:,3), sX(:,:,4), sX(:,:,5), sX(:,:,6), sX(:,:,7), sX(:,:,8), sX(:,:,9)] = ...
    arrayfun(@stencil, x(1:end-8,:), x(2:end-7,:), x(3:end-6,:), x(4:end-5,:), ...
    x(5:end-4,:), x(6:end-3,:), x(7:end-2,:), x(8:end-1,:), x(9:end,:), methodsX);

[sY(:,:,1), sY(:,:,2), sY(:,:,3), sY(:,:,4), sY(:,:,5), sY(:,:,6), sY(:,:,7), sY(:,:,8), sY(:,:,9)] = ...
    arrayfun(@stencil, y(1:end-8,:), y(2:end-7,:), y(3:end-6,:), y(4:end-5,:), ...
    y(5:end-4,:), y(6:end-3,:), y(7:end-2,:), y(8:end-1,:), y(9:end,:), methodsY);

sY = -permute(sY, [2,1,3]);
sX = -sX;

sY(isnan(sY) | isinf(sY)) = 0;
sX(isnan(sX) | isinf(sX)) = 0;

s_X.(x_dims{1}) = sX;
s_X.(x_dims{2}) = sY;

if nargout == 2
    nstencil = size(sX);
    nstencil = nstencil(end);
end
end

function [s1, s2, s3, s4, s5, s6, s7, s8, s9] = stencil(x1, x2, x3, x4, x5, x6, x7, x8, x9, method)
switch method
    case 0
        % In boundary
        s1 = 0; s2 = 0; s3 = 0; s4 = 0; s5 = 0; s6 = 0; s7 = 0; s8 = 0; s9 = 0;
    case 1 
        % Forward 1st order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = firstOrder(x5, x6, 1);
    case 2 
        % Forward 2nd order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = secondOrder(x5, x6, x7, 1);
    case 3
        % Forward 3rd order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = thirdOrder(x5, x6, x7, x8, 1);
    case 4
        % Forward 4th order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = fourthOrder(x5, x6, x7, x8, x9, 1);
    case 5
        % Forward biased 3rd order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = thirdOrder(x4, x5, x6, x7, 2);
    case 6
        % Forward biased 4th order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = fourthOrder(x4, x5, x6, x7, x8, 2);
    case 7
        % 1 point stencil
        s1 = 0; s2 = 0; s3 = 0; s4 = 0; s5 = 0; s6 = 0; s7 = 0; s8 = 0; s9 = 0;
    case 8
        % Central 2nd order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = secondOrder(x4, x5, x6, 2);
    case 9
        % Central 4th order
        [s1, s2, s3, s4, s5, s6, s7, s8, s9] = fourthOrder(x3, x4, x5, x6, x7, 3);
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