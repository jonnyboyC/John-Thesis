function [gmm_sim, km_sim] = classify_model(km, gmm, modal_amp_model, num_clusters, outlier_mode)
% CLASSIFY_MODEL use the scheme described in John Chabot's Thesis 2015 to
% classify each point in a simulated model to those of the orignal data
%
%   [gmm, km] = CLASSIFY_MODEL(km, gmm, modal_amp_model)

[km_sim.groups, distance] = knnsearch(km.centers, modal_amp_model);
[gmm_sim.groups, ~, ~, ~, mahal_distance] = cluster(gmm.models, modal_amp_model);

if ~isempty(outlier_mode.km) 
    % Location of outlier mode
    outlier = num_clusters + 1;
    
    % outlier modes in k-means
    center_dis = outlier_mode.km*max(max(pdist2(km.centers, km.centers)));
    km_sim.groups(distance > center_dis) = outlier;
    
    % outlier modes in gmm
    gauss_sigma = outlier_mode.gmm;
    mahal_distance = min(mahal_distance,[],2);
    gmm_sim.groups(mahal_distance > gauss_sigma) = outlier;
end
end