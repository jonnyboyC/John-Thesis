function [x, y, u, v, num_processed] = load_jet_dat(num_images, direct)
% LOAD_JET_DAT load required dat files to run the DNS simulated axysmetric
% data
%
%   [x, y, u, v] = LOAD_JET_DAT(img_files, num_files, num_images,
%   image_range, flip, direct, see help POD_GEN for informaiton

[x, y, ~, num_zones] = load_grid_dat(direct);
[u, v, ~] = load_velocity_dat(direct, num_images, num_zones);
num_processed = num_images;
end