function [gmm, km] = gen_clusters(modal_amp, modes, num_clusters, num_cores)
% If cluster has not already been generated produce cluster

valid = false;
multiplier = 1;

% produce clusters
[km.groups, km.centers, gmm.models, gmm.groups] ...
    = cluster_amp(modal_amp, modes, num_clusters, num_cores);

% get secondary information
km.stoch  = gen_stochastic_matrix(num_clusters, km.groups, multiplier, valid);
gmm.stoch = gen_stochastic_matrix(num_clusters, gmm.groups, multiplier, valid);

gmm.prob = calc_probability(gmm.stoch, gmm.groups);
km.prob  = calc_probability(km.stoch, km.groups);
end