function coefficient = surf_inner_prod(var1, var2, volume, bound)
num1 = size(var1, 2);
num2 = size(var2, 2);

bound(bound > 0) = 1;
bound = logical(bound);
surf_vol = volume(bound);

mask1 = repmat(reshape(bound, numel(bound), 1), 1, num1);
mask2 = repmat(reshape(bound, numel(bound), 1), 1, num2);

surf1 = var1(mask1);
surf2 = var2(mask2);

surf1 = reshape(surf1, [], num1);
surf2 = reshape(surf2, [], num2);

coefficient = (surf1'*(surf2.*(surf_vol*ones(1,num2))))';
end

