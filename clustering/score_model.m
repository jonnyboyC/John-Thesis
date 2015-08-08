function [frob_km, frob_gm, prob_km, prob_gm, completed] = score_model(score_info)
%  SCORE_GAL generate model scores based on method outline in John Chabot's
%  Thesis 2015 classify simulated data then apply score metrics. If
%  clusters were not precompiled from POD_Gen will generate them here
%
% [frob_km, frob_gm, prob_km, prob_gm, completed] = SCORE_GAL(km, gmm, ...
%   integration, tspan, num_clusters, mutiplier, int_time, it, direct)

% unpack variables
km              = score_info.km;
gmm             = score_info.gmm;
integration     = score_info.integration;
tspan           = score_info.tspan;
modal_amp       = score_info.modal_amp;
num_clusters    = score_info.num_clusters;
int_time        = score_info.int_time;
num_cores       = score_info.num_cores;
modes           = score_info.modes;
custom          = score_info.custom;
direct          = score_info.direct;
multiplier      = score_info.multiplier;
MOD             = score_info.MOD;

% Free memory
clear score_info

% inputs for gen_stochoastic_matrix
save_figures = {};

% Get eddy visocity models
m = flow_comps(integration.modal_amp);
m_comps = flow_ncomps(integration.modal_amp);

h_km = figure;
h_gm = figure;

% Strip off mean mode if we need to cluster
if ~MOD
   modal_amp = modal_amp(:,2:end); 
end

% Select or generate clusters from data
if custom
    [gmm, km] = gen_clusters(modal_amp, modes, num_clusters, num_cores);
else
    if (length(modes) - 1) <= length(km)
        idx = length(modes) - 1;
        km = km{idx};
        gmm = gmm{idx};
    else
       [gmm, km] = gen_clusters(modal_amp, modes, num_clusters, num_cores); 
    end
end
    

valid_transition = true;

% Loop through each model calculating scores
for j = 1:m_comps
    
    s = flow_comps(integration.modal_amp.(m{j}));
    s_comps = flow_ncomps(integration.modal_amp.(m{j}));
    
    for k = 1:s_comps
        
        % If something happened the integration don't attempt to classify
        if isempty(integration.t.(m{j}).(s{k})) || isfield(integration.t.(m{j}).(s{k}), 'error')
            frob_km.(m{j}).(s{k}) = [];
            frob_gm.(m{j}).(s{k}) = [];
            
            prob_km.(m{j}).(s{k}) = [];
            prob_gm.(m{j}).(s{k}) = [];
            
            completed.(m{j}).(s{k}) = false;
            
            fprintf('%s %s took longer than %d seconds to integrate and was discarded\n', ...
                strrep((m{j}), '_', ' '), strrep((s{k}), '_', ' '), int_time);
            continue;
        end
        
        % Classify each point to a cluster
        if MOD
            km_sim.groups = knnsearch(km.centers, integration.modal_amp.(m{j}).(s{k}));
            gmm_sim.groups = cluster(gmm.models, integration.modal_amp.(m{j}).(s{k}));
        else
            km_sim.groups = knnsearch(km.centers, integration.modal_amp.(m{j}).(s{k})(:,2:end));
            gmm_sim.groups = cluster(gmm.models, integration.modal_amp.(m{j}).(s{k})(:,2:end));
        end
        
        % Generate a stochastic matrix from transitions
        km_sim.stoch = gen_stochastic_matrix(num_clusters, km_sim.groups, multiplier, valid_transition);
        gmm_sim.stoch = gen_stochastic_matrix(num_clusters, gmm_sim.groups, multiplier, valid_transition);
        
        % log probability of observing the simulate chain using MLE model
        km_sim.prob = calc_probability(km_sim.stoch, km.groups);
        gmm_sim.prob = calc_probability(gmm_sim.stoch, gmm.groups);
        
        % Frobenius norm of transition matrices
        frob_km.(m{j}).(s{k}) = norm(km_sim.stoch - km.stoch, 'fro');
        frob_gm.(m{j}).(s{k}) = norm(gmm_sim.stoch - gmm.stoch, 'fro');
        
        % Relative likelihood between observed chains
        prob_km.(m{j}).(s{k}) = km_sim.prob - km.prob;
        prob_gm.(m{j}).(s{k}) = gmm_sim.prob - gmm.prob;
        
        plot_stochastic_matrix(km_sim, save_figures, direct, h_km);
        plot_stochastic_matrix(gmm_sim, save_figures, direct, h_gm);
        
        % Temporary reporting
        if integration.t.(m{j}).(s{k})(end) == tspan(end)
            fprintf('%s %s finished got km_frob = %2.4f and gm_frob = %2.4f and km_prob = %2.4f and gm_prob = %2.4f\n', ...
                strrep((m{j}), '_', ' '), strrep((s{k}), '_', ' '),frob_km.(m{j}).(s{k}), ...
                frob_gm.(m{j}).(s{k}), prob_km.(m{j}).(s{k}), prob_gm.(m{j}).(s{k}));
            completed.(m{j}).(s{k}) = true;
        else
            fprintf('%s %s diverged got km_frob = %2.4f and gm_frob = %2.4f and km_prob = %2.4f and gm_prob = %2.4f\n', ...
                strrep((m{j}), '_', ' '), strrep((s{k}), '_', ' '),frob_km.(m{j}).(s{k}), ...
                frob_gm.(m{j}).(s{k}), prob_km.(m{j}).(s{k}), prob_gm.(m{j}).(s{k}));
            completed.(m{j}).(s{k}) = false;
        end
    end
end
end

