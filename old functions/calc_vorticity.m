function vorticity = calc_vorticity(data, u, v, dimensions, bnd_idx)
%% TODO do this more efficiently 

% Figure out how many loops are needed
number2calc = size(u, 2);

% Precalculate values for finite difference
[xxi, yxi, xet, yet, aj] = metric2(data.xg, data.yg);
z = ones(dimensions);

% Calculate coefficients for u's & v's derivatives
udy = zeros(dimensions(1), dimensions(2), number2calc);
vdx = zeros(dimensions(1), dimensions(2), number2calc);

parfor i = 1:number2calc
    
    [dxic_x, detc_x] = visder2(reshape(u(:,i), size(z,1), size(z,2)),...
        dimensions(1), dimensions(2), z, bnd_idx);
     
    [dxic_y, detc_y] = visder2(reshape(v(:,i), size(z,1), size(z,2)),...
        dimensions(1), dimensions(2), z, bnd_idx);
    
    udy(:,:,i) = aj.*(detc_x.*xxi-dxic_x.*xet);
    vdx(:,:,i) = aj.*(dxic_y.*yet-detc_y.*yxi);
end

% Calculate vorticity
vorticity = vdx - udy;

vorticity = reshape(vorticity, numel(data.xg), size(u,2));
end