function coefficient = surf_inner_prod(var1, var2, vol_frac, bound)
mask = repmat(reshape(bound, numel(bound), 1), 1, size(var1, 2));
surf1 = var1(mask);
mask = repmat(reshape(bound, numel(bound), 1), 1, size(var2, 2));
surf2 = var2(mask);
coefficient = 0;
end

