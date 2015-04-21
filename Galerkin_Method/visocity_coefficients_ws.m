function l = visocity_coefficients_ws(coef_problem)
% Unpack Variables
x           = coef_problem.x;
y           = coef_problem.y;
pod_u       = coef_problem.pod_u;
pod_v       = coef_problem.pod_v;
dimensions  = coef_problem.dimensions;
vol_frac    = coef_problem.vol_frac;
run_num     = coef_problem.run_num;
direct      = coef_problem.direct;
override_coef = coef_problem.override_coef;
bnd_idx     = coef_problem.bnd_idx;
bnd_x       = coef_problem.bnd_x;
bnd_y       = coef_problem.bnd_y;

clear coef_problem

num_elem = numel(x);
num_modes = size(pod_u, 2);
bnd_x = reshape(bnd_x, num_elem, 1);
bnd_y = reshape(bnd_y, num_elem, 1);
[bnd_x, bnd_y] = strip_boundaries(bnd_idx, bnd_x, bnd_y);

% If we are using the same number of cutoff modes and overwrite is set to
% false look for previous data
if override_coef == false;
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
[pod_udx, pod_udy, pod_vdx, pod_vdy, pod_u, pod_v, vol_frac] = ...
    components_ws(x, y, pod_u, pod_v, dimensions, vol_frac, num_modes, bnd_idx, bnd_x, bnd_y);


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
save([direct '\Viscous Coeff\Coeff_' num2str(run_num) '_wk_m' num2str(num_modes) '.mat'], ...
    'l', 'cutoff', 'run_num', '-v7.3'); 

end
