function model = stationary(model)
% STATIONARY Compute the stationary distribution of an ergodic Markov chain
% caracterized by its stochastic matrix of cluster model such that
% P = Pr(x_{k} | x_{k - 1}) , sum(P , 2) = 1. Add on additional field to
% model.
%
% model = stationary(model);
%

if ~isfield(model, 'stoch')
    error('Provide a k-mean or gmm model see gen_data_clusters');
end

P = model.stoch;

[M, N] = size(P);

if (M ~= N)
    error('Model contains an invalid stochastic matrix')
end

% Stationary distrubtion is eigenvalue decomposition of stochastic matrix
[V , D] = eig(P');
[~, idx] = sort(diag(D));
 
% Normalize
model.stat = (V(:, idx(end))/sum(V(:, idx(end))))';