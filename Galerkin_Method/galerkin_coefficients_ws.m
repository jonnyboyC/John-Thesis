function l = galerkin_coefficients_ws(coef_problem)
% GALERKIN_COEFFICIENT_WS calculate the linear terms resultant from the
% visocity term in the Navier Stoke equations using Green's Identiy
%
% l = GALERKIN_COEFFICIENT_WS(coef_problem)

% Unpack Variables
X           = coef_problem.X;
pod_U       = coef_problem.pod_U;
dimensions  = coef_problem.dimensions;
volume      = coef_problem.volume;
run_num     = coef_problem.run_num;
direct      = coef_problem.direct;
override_coef = coef_problem.override_coef;
bnd_idx     = coef_problem.bnd_idx;
bnd_X       = coef_problem.bnd_X;
custom      = coef_problem.custom;
uniform     = coef_problem.uniform;


clear coef_problem

[x, u] = flow_comps_ip(X, pod_U);
dims = flow_dims(X);

num_elem = numel(X.(x{1}));
num_modes = size(pod_U.(u{1}), 2);

for i = 1:dims
   bnd_X.(x{i}) = reshape(bnd_X.(x{i}), num_elem, 1); 
end

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
[pod_UdX, pod_U, bnd_X, volume] = ...
    components_ws(X, pod_U, dimensions, volume, num_modes, num_elem, uniform, bnd_idx, bnd_X);

% Replace Lapcalian term with weak formulation by green's identity 
l_weak_volume = 0;
l_weak_surf = 0;
for i = 1:dims
    for j = 1:dims
        l_weak_volume = l_weak_volume ...
            - inner_prod(pod_UdX.(u{i}).(x{j}), pod_UdX.(u{i}).(x{j}), volume);
        l_weak_surf = l_weak_surf ...
            + surf_inner_prod(pod_UdX.(u{i}).(x{j}), pod_U.(u{i}), volume, bnd_X.(x{j}));
    end
end

l = l_weak_volume + l_weak_surf;

cutoff = num_modes;
if ~custom
    save([direct '\Viscous Coeff\Coeff_' num2str(run_num) '_wk_m' num2str(num_modes) '.mat'], ...
        'l', 'cutoff', 'run_num', '-v7.3');
end

end
