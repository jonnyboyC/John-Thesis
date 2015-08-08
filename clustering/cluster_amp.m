function [km_groups, km_centers, gmm_models, gmm_groups] = cluster_amp(modal_amp, modes, num_clusters, varargin)
% CLUSTER_AMP produces clusters for k-mean and gaussian mixture model for k
% clusters
%
%   [km_groups, km_centers, gm_models, gm_groups] ...
%       = CLUSTER_AMP(modal_amp, modes, num_clusters)
%
%   [km_groups, km_centers, gm_models, gm_groups] ...
%       = CLUSTER_AMP(modal_amp, modes, num_clusters, num_cores)

if nargin == 3
    num_cores = 2;
elseif nargin == 4
    num_cores = varargin{1};
else
    error('Only 4 inputs, type help cluster_amp for details');
end

% Clustering options
if num_cores ~= 1
    options = statset('UseParallel', true);
else
    options = statset('UseParallel', false);
end
gm_options = statset('MaxIter', 1000);

% clustering random start replicates
k_replicate = 50;
gm_replicate = 10;

% Cluster based on k-mean
[km_groups, km_centers] = kmeans(modal_amp(:,modes), ...
    num_clusters, 'Replicates', k_replicate, 'Options', options);

% Fit gaussian mixture model
try 
    % First attempt to produce a fit with independent covariance matrix
    gmm_models = fitgmdist(modal_amp(:,modes), num_clusters, ...
        'Replicates', gm_replicate, 'SharedCov', false, 'Options', gm_options);
catch error_msg
    % If it fails fall backed to shared covariance matrices
    if (strcmp(error_msg.identifier, 'stats:gmdistribution:IllCondCovAllReps'))
        gmm_models = fitgmdist(modal_amp(:,modes), num_clusters, ...
            'Replicates', gm_replicate, 'SharedCov', true, 'Options', gm_options);
    else
        rethrow(error_msg);
    end
end

% Cluster based on gaussian mixture model
gmm_groups = cluster(gmm_models, modal_amp(:,modes));
end