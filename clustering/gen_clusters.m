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
[km, gmm] = cluster_amp(modal_amp, modes, num_clusters, num_cores);

% get secondary information
km  = gen_stochastic_matrix(km, num_clusters, multiplier, classify, outlier_mode);
gmm = gen_stochastic_matrix(gmm, num_clusters,  multiplier, classify, outlier_mode);

% Calculate exact likelihood of observing Markov Chain for the MLE
km = calc_likelihood(km);
gmm  = calc_likelihood(gmm);
end