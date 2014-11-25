function vorticity = calc_vorticity(data, u, v, dimensions, bnd_idx)
%% Calculate vorticity in the flow
[xxi, yxi, xet, yet, aj] = metric2(data.xg, data.yg);
z = ones(dimensions);

% Calculate coefficients for u's & v's derivatives
[~, ~, udy, ~] = derivatives(u, dimensions, z, xxi, yxi,...
    xet, yet, aj, bnd_idx);
[vdx, ~, ~, ~] = derivatives(v, dimensions, z, xxi, yxi,...
    xet, yet, aj, bnd_idx);

% Calculate vorticity
vorticity = vdx - udy;
vorticity = reshape(vorticity, numel(data.xg), size(u,2));
end