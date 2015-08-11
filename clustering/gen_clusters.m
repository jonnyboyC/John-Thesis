function [gmm, km] = gen_clusters(modal_amp, modes, num_clusters, num_cores, outlier_mode)
% GEN_CLUSTERS generates clusters from empirical or numerical data using
% k-means and gaussian mixture models. Additionally calculates the estimate
% of the stochastic matrix and the probability of obsvering the data Markov
% chain given the estimated stochastic matrix
%
%   [gmm, km] = GEN_CLUSTERS(modal_amp, modes, num_clusters, num_cores)

classify = false;
multiplier = 1;

% produce clusters
[km.groups, km.centers, gmm.models, gmm.groups] ...
    = cluster_amp(modal_amp, modes, num_clusters, num_cores);

% get secondary information
km.stoch  = gen_stochastic_matrix(num_clusters, km.groups, multiplier, classify, outlier_mode);
gmm.stoch = gen_stochastic_matrix(num_clusters, gmm.groups, multiplier, classify, outlier_mode);

% Calculate exact likelihood of observing Markov Chain for the MLE
gmm.like = calc_likelihood(gmm.stoch, gmm.groups);
km.like  = calc_likelihood(km.stoch, km.groups);
end