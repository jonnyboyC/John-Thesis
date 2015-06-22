function l = galerkin_coefficients_ws(coef_problem)
% GALERKIN_COEFFICIENT_WS calculate the linear terms resultant from the
% visocity term in the Navier Stoke equations using Green's Identiy
%
% l = GALERKIN_COEFFICIENT_WS(coef_problem)

% Unpack Variables
X           = coef_problem.X;
pod_U       = coef_problem.pod_U;
dimensions  = coef_problem.dimensions;
vol_frac    = coef_problem.vol_frac;
run_num     = coef_problem.run_num;
direct      = coef_problem.direct;
override_coef = coef_problem.override_coef;
bnd_idx     = coef_problem.bnd_idx;
bnd_X       = coef_problem.bnd_X;
custom      = coef_problem.custom;

clear coef_problem

[x, u] = flow_comps_ip(X, pod_U);
dims = flow_dims(X);

num_elem = numel(X.(x{1}));
num_modes = size(pod_U.(u{1}), 2);

[bnd_x, bnd_y] = strip_boundaries(bnd_idx, bnd_x, bnd_y);

% If we are using the same number of cutoff modes and overwrite is set to
% false look for previous data
if override_coef == false && custom == false
   saved_files = dir([direct '\Viscous Coeff\Coeff_*']);
   if size(saved_files,1) ~= 0
       match_run = regexp({saved_files.name}, num2str(run_num));
       match_modes = regexp({saved_files.name}, ['m' num2str(num_modes) '\.']);
       weak = regexp({saved_files.name}, 'wk');
       if any(~cellfun(@isempty, match_run) & ~cellfun(@isempty, match_modes) & ~cellfun(@isempty, weak))
           data = load([direct '\Viscous Coeff\Coeff_' num2str(run_num) '_wk_m' num2str(num_modes) '.mat']);
           l = data.l;
           return;
       end
   end
end

% Calculate terms, allows for nonuniform mesh
[pod_UdX, pod_U, vol_frac] = ...
    components_ws(X, pod_U, dimensions, vol_frac, num_modes, num_elem, bnd_idx, bnd_X);


% Calculated Weak Solution Gradient inner products
ccux = -inner_prod(pod_udx, pod_udx, vol_frac);
ccvx = -inner_prod(pod_vdx, pod_vdx, vol_frac);

ccuy = -inner_prod(pod_udy, pod_udy, vol_frac);
ccvy = -inner_prod(pod_vdy, pod_vdy, vol_frac);

% Calulate Weak Solution Surface Integral inner products
surf_cu = surf_inner_prod(pod_udx, pod_u, vol_frac, bnd_x) + ...
          surf_inner_prod(pod_udy, pod_u, vol_frac, bnd_y);
surf_cv = surf_inner_prod(pod_vdx, pod_v, vol_frac, bnd_x) + ...
          surf_inner_prod(pod_vdy, pod_v, vol_frac, bnd_y);

l = ccux + ccuy + ccvx + ccvy + surf_cu + surf_cv;

cutoff = num_modes;
if ~custom
    save([direct '\Viscous Coeff\Coeff_' num2str(run_num) '_wk_m' num2str(num_modes) '.mat'], ...
        'l', 'cutoff', 'run_num', '-v7.3');
end

end
