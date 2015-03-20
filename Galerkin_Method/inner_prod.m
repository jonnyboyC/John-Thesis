function coefficient = inner_prod(var1, var2, vol_frac)
% Take the inner product, producing a portion of the energy term
coefficient = (var1'*(var2.*(repmat(vol_frac,1, size(var2,2)))))';
end