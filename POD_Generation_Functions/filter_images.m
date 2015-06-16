function [U, mean_U] = filter_images(U, bnd_idx, num_images, ensemble_dim)
% FILTER_IMAGES apply guided filter to the flow images, currently is not
% ideal because of some smooth at the boundary, could use some rework
%
%   [U, mean_U] = FILTER_IMAGES(U, bnd_idx, num_images, ensemble_dim)

% Get flow information
comps = flow_ncomps(U);
dims = flow_dims(U);
u = flow_comps(U);

% Generate index
idx = flow_index({1 1}, dims(end), U);

% Perform filtering
for i = 1:num_images
    for j = 1:comps
        idx{end} = i;
        U_temp = U.(u{j})(idx{:});
        U_gauss= imgaussfilt(U_temp);
        U_gauss = imgaussfilt(U_gauss);
        U_gauss(bnd_idx == 1) = U_temp(bnd_idx == 1);
        U.(u{j})(idx{:}) = imguidedfilter(U_gauss);
    end
end

% Update mean values
mean_U = mean_comps(U, ensemble_dim);