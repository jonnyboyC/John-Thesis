function [gm, km] = gen_data_clusters(modal_amp, modes, num_clusters, num_cores)
% GEN_DATA_CLUSTER generate kmean and gmm cluster models from an empirical data
% source
%
% [gm, km] = GEN_DATA_CLUSTERS(modal_amp, modes, num_clusters, num_cores)

valid = false;
multiplier = 1;

% produce clusters
[km.groups, km.centers, gm.models, gm.groups] ...
    = cluster_amp(modal_amp, modes, num_clusters, num_cores);

% Generate stochastic matrices for both clustering methods
km.stoch = gen_stochastic_matrix(num_clusters, km.groups, multiplier, valid);
gm.stoch = gen_stochastic_matrix(num_clusters, gm.groups, multiplier, valid);

% Calculate probability of observing data's chain
gm.prob = calc_probability(gm.stoch, gm.groups);
km.prob = calc_probability(km.stoch, km.groups);
end