function [bnd_idx] = flow_boundaries(mean_u, mean_v)
% Make all none zero mean velocity 1
bnd_idx = double((mean_u + mean_v ~= 0));

% Intially make points with zero mean velocity -1
bnd_idx(bnd_idx == 0) = -1;

% Use built in edge detection to boundary points change points to 0
bnd_idx(edge(bnd_idx, 'canny')) = 0;
end