function coefficient = surf_inner_prod(var1, var2, volume, bound)
% SURF_INNER_PROD take the inner product for the surface intergral when
% using the weak formulation of Navier-Stokes
%
%   coefficients = SURF_INNER_PROD(var1, var2, volume, bound)

% TODO currently take inner produce for a 1 pixel volume along the
% bndaries, could probably be made more accuracy by somehow calculating the
% actual line/surface volume

% Determine number of images
images1 = size(var1, 2);
images2 = size(var2, 2);

% Find values in boundaries
bound(bound > 0) = 1;
bound = logical(bound);
surf_vol = volume(bound);

% Mask and strip points not in mask
mask1 = repmat(reshape(bound, numel(bound), 1), 1, images1);
mask2 = repmat(reshape(bound, numel(bound), 1), 1, images2);

surf1 = var1(mask1);
surf2 = var2(mask2);

surf1 = reshape(surf1, [], images1);
surf2 = reshape(surf2, [], images2);

% Take inner produce of surface
coefficient = (surf1'*(surf2.*(surf_vol*ones(1,images2))))';
end

