function relations = cluster_plots(clusters_info)
% Will be garbage code :(

results_scores  = clusters_info.results_scores;
TKE             = clusters_info.TKE;
frequency       = clusters_info.frequency;
phase_comp      = clusters_info.phase_comp;
amplitude1      = clusters_info.amplitude1;
amplitude2      = clusters_info.amplitude2;
num_clusters    = clusters_info.num_clusters;
direct          = clusters_info.direct;

% Field names
tke_fields  = {'mean', 'std', 'median'};
amp1_fields = {'mean', 'std', 'median'};
amp2_fields = {'mean', 'std', 'median'};
freq_fields = {'amp1', 'amp2' 'amp3'};
phase_fields= {'mean', 'std', 'median'};

relations = struct;

[fkm, full_name]    = merge_struct(results_scores,  'frob_km');
fgmm                = merge_struct(results_scores,  'frob_gmm');
pkm                 = merge_struct(results_scores,  'like_km');
pgmm                = merge_struct(results_scores,  'like_gmm');

% If literally no models stayed bounded
if isempty(full_name)
    return;
end

% Get continous lists of models and values
tke.mean        = merge_struct(TKE, 'mean_diff');
tke.std         = merge_struct(TKE, 'std_diff');
tke.median      = merge_struct(TKE, 'median_diff');

amp1.mean       = merge_struct(amplitude1, 'mean_diff');
amp1.std        = merge_struct(amplitude1, 'std_diff');
amp1.median     = merge_struct(amplitude1, 'median_diff');

amp2.mean       = merge_struct(amplitude2, 'mean_diff');
amp2.std        = merge_struct(amplitude2, 'std_diff');
amp2.median     = merge_struct(amplitude2, 'median_diff');

if ~isempty(frequency)
    fft.amp1    = merge_struct(frequency, 'amp1_diff');
    fft.amp2    = merge_struct(frequency, 'amp2_diff');
    fft.amp3    = merge_struct(frequency, 'amp3_diff');
end

if ~isempty(phase_comp)
    phase.mean  = merge_struct(phase_comp, 'mean_diff');
    phase.std   = merge_struct(phase_comp, 'std_diff');
    phase.median= merge_struct(phase_comp, 'median_diff');
end

sub_name = cell(size(full_name));
for i = 1:length(full_name)
    underscore = strfind(full_name{i}, '_');
    sub_name{i} = full_name{i}(underscore+1:end);
end

%% TKE relations

for i = 1:length(tke_fields)
    
    % compare to k-means trasnition norm score
    [~, R2, R, p] = score_plot(full_name, sub_name, fkm, tke.(tke_fields{i}), 'k-means \sigma_{d}', ...
        ['TKE ' tke_fields{i}], ['TKE ' tke_fields{i} ' vs. k-mean \sigma_{d}'], direct, ...
        num_clusters, ['fkm_v_tke_' tke_fields{i}], 'log');
    
	relations.fkm.tke.(tke_fields{i}).R2 = R2;
 	relations.fkm.tke.(tke_fields{i}).R = R;
    relations.fkm.tke.(tke_fields{i}).p = p;
    
    % compare to gmm trasnition norm score
    [~, R2, R, p] = score_plot(full_name, sub_name, fgmm, tke.(tke_fields{i}), 'GMM \sigma_{d}', ...
        ['TKE ' tke_fields{i}], ['TKE ' tke_fields{i} ' vs. GMM \sigma_{d}'], direct, ...
        num_clusters, ['fgmm_v_tke_' tke_fields{i}], 'log');
    
	relations.fgmm.tke.(tke_fields{i}).R2 = R2;
 	relations.fgmm.tke.(tke_fields{i}).R = R;
    relations.fgmm.tke.(tke_fields{i}).p = p;
    
    % compare to k-mean likelihood score
    [~, R2, R, p] = score_plot(full_name, sub_name, pkm, tke.(tke_fields{i}), 'k-means \sigma_{l}', ...
        ['TKE ' tke_fields{i}], ['TKE ' tke_fields{i} ' vs. k-mean \sigma_{l}'], direct, ...
        num_clusters, ['pkm_v_tke_' tke_fields{i}], 'log');
    
	relations.pkm.tke.(tke_fields{i}).R2 = R2;
 	relations.pkm.tke.(tke_fields{i}).R = R;
    relations.pkm.tke.(tke_fields{i}).p = p;
    
    % compare to gmm likelihood score
    [~, R2, R, p] = score_plot(full_name, sub_name, pgmm, tke.(tke_fields{i}), 'GMM \sigma_{l}', ...
        ['TKE ' tke_fields{i}], ['TKE ' tke_fields{i} ' vs. GMM \sigma_{l}'], direct, ...
        num_clusters, ['pgmm_v_tke_' tke_fields{i}], 'log');
    
	relations.pgmm.tke.(tke_fields{i}).R2 = R2;
 	relations.pgmm.tke.(tke_fields{i}).R = R;
    relations.pgmm.tke.(tke_fields{i}).p = p;
end

%% Amplitude 1 relations

for i = 1:length(amp1_fields)
    
    % compare to k-means trasnition norm score
    [~, R2, R, p] = score_plot(full_name, sub_name, fkm, amp1.(amp1_fields{i}), 'k-means \sigma_{d}', ...
        ['Amplitude a_{1} ' amp1_fields{i}], ['Amplitude a_{1} ' amp1_fields{i} ' vs. k-mean \sigma_{d}'], direct, ...
        num_clusters, ['fkm_v_amp1_' amp1_fields{i}], 'log');
    
	relations.fkm.amp1.(amp1_fields{i}).R2 = R2;
 	relations.fkm.amp1.(amp1_fields{i}).R = R;
    relations.fkm.amp1.(amp1_fields{i}).p = p;
    
    % compare to gmm trasnition norm score
    [~, R2, R, p] = score_plot(full_name, sub_name, fgmm, amp1.(amp1_fields{i}), 'GMM \sigma_{d}', ...
        ['Amplitude a_{1} ' amp1_fields{i}], ['Amplitude a_{1} ' amp1_fields{i} ' vs. GMM \sigma_{d}'], direct, ...
        num_clusters, ['fgmm_v_amp1_' amp1_fields{i}], 'log');
    
	relations.fgmm.amp1.(amp1_fields{i}).R2 = R2;
 	relations.fgmm.amp1.(amp1_fields{i}).R = R;
    relations.fgmm.amp1.(amp1_fields{i}).p = p;
    
    % compare to k-mean likelihood score
    [~, R2, R, p] = score_plot(full_name, sub_name, pkm, amp1.(amp1_fields{i}), 'k-means \sigma_{l}', ...
        ['Amplitude a_{1} ' amp1_fields{i}], ['Amplitude a_{1} ' amp1_fields{i} ' vs. k-mean \sigma_{l}'], direct, ...
        num_clusters, ['pkm_v_amp1_' amp1_fields{i}], 'log');
    
	relations.pkm.amp1.(amp1_fields{i}).R2 = R2;
 	relations.pkm.amp1.(amp1_fields{i}).R = R;
    relations.pkm.amp1.(amp1_fields{i}).p = p;
    
    % compare to gmm likelihood score
    [~, R2, R, p] = score_plot(full_name, sub_name, pgmm, amp1.(amp1_fields{i}), 'GMM \sigma_{l}', ...
        ['Amplitude a_{1} ' amp1_fields{i}], ['Amplitude a_{1} ' amp1_fields{i} ' vs. GMM \sigma_{l}'], direct, ...
        num_clusters, ['pgmm_v_amp1_' amp1_fields{i}], 'log');
    
	relations.pgmm.amp1.(amp1_fields{i}).R2 = R2;
 	relations.pgmm.amp1.(amp1_fields{i}).R = R;
    relations.pgmm.amp1.(amp1_fields{i}).p = p;
end

%% Amplitude 2 relations

for i = 1:length(amp2_fields)
    
    % compare to k-means trasnition norm score
    [~, R2, R, p] = score_plot(full_name, sub_name, fkm, amp2.(amp2_fields{i}), 'k-means \sigma_{d}', ...
        ['Amplitude a_{2} ' amp2_fields{i}], ['Amplitude a_{2} ' amp2_fields{i} ' vs. k-mean \sigma_{d}'], direct, ...
        num_clusters, ['fkm_v_amp2_' amp2_fields{i}], 'log');
    
	relations.fkm.amp2.(amp2_fields{i}).R2 = R2;
 	relations.fkm.amp2.(amp2_fields{i}).R = R;
    relations.fkm.amp2.(amp2_fields{i}).p = p;
    
    % compare to gmm trasnition norm score
    [~, R2, R, p] = score_plot(full_name, sub_name, fgmm, amp2.(amp2_fields{i}), 'GMM \sigma_{d}', ...
        ['Amplitude a_{2} ' amp2_fields{i}], ['Amplitude a_{2} ' amp2_fields{i} ' vs. GMM \sigma_{d}'], direct, ...
        num_clusters, ['fgmm_v_amp2_' amp2_fields{i}], 'log');
    
	relations.fgmm.amp2.(amp2_fields{i}).R2 = R2;
 	relations.fgmm.amp2.(amp2_fields{i}).R = R;
    relations.fgmm.amp2.(amp2_fields{i}).p = p;
    
    % compare to k-mean likelihood score
    [~, R2, R, p] = score_plot(full_name, sub_name, pkm, amp2.(amp2_fields{i}), 'k-means \sigma_{l}', ...
        ['Amplitude a_{2} ' amp2_fields{i}], ['Amplitude a_{2} ' amp2_fields{i} ' vs. k-mean \sigma_{l}'], direct, ...
        num_clusters, ['pkm_v_amp2_' amp2_fields{i}], 'log');
    
	relations.pkm.amp2.(amp2_fields{i}).R2 = R2;
 	relations.pkm.amp2.(amp2_fields{i}).R = R;
    relations.pkm.amp2.(amp2_fields{i}).p = p;
    
    % compare to gmm likelihood score
    [~, R2, R, p] = score_plot(full_name, sub_name, pgmm, amp2.(amp2_fields{i}), 'GMM \sigma_{l}', ...
        ['Amplitude a_{2} ' amp2_fields{i}], ['Amplitude a_{2} ' amp2_fields{i} ' vs. GMM \sigma_{l}'], direct, ...
        num_clusters, ['pgmm_v_amp2_' amp2_fields{i}], 'log');
    
	relations.pgmm.amp2.(amp2_fields{i}).R2 = R2;
 	relations.pgmm.amp2.(amp2_fields{i}).R = R;
    relations.pgmm.amp2.(amp2_fields{i}).p = p;
end

%% Frequency Relations

if ~isempty(frequency)
for i = 1:length(freq_fields)
    
    % compare to k-means trasnition norm score
    [~, R2, R, p] = score_plot(full_name, sub_name, fkm, fft.(freq_fields{i}), 'k-means \sigma_{d}', ...
        ['Frequency Offset ' freq_fields{i}], ['Frequency Offset ' freq_fields{i} ' vs. k-mean \sigma_{d}'], direct, ...
        num_clusters, ['fkm_v_freq_' freq_fields{i}], 'log');
    
	relations.fkm.fft.(freq_fields{i}).R2 = R2;
 	relations.fkm.fft.(freq_fields{i}).R = R;
    relations.fkm.fft.(freq_fields{i}).p = p;
    
    % compare to gmm trasnition norm score
    [~, R2, R, p] = score_plot(full_name, sub_name, fgmm, fft.(freq_fields{i}), 'GMM \sigma_{d}', ...
        ['Frequency Offset ' freq_fields{i}], ['Frequency Offset ' freq_fields{i} ' vs. GMM \sigma_{d}'], direct, ...
        num_clusters, ['fgmm_v_freq_' freq_fields{i}], 'log');
    
	relations.fgmm.fft.(freq_fields{i}).R2 = R2;
 	relations.fgmm.fft.(freq_fields{i}).R = R;
    relations.fgmm.fft.(freq_fields{i}).p = p;
    
    % compare to k-mean likelihood score
    [~, R2, R, p] = score_plot(full_name, sub_name, pkm, fft.(freq_fields{i}), 'k-means \sigma_{l}', ...
        ['Frequency Offset ' freq_fields{i}], ['Frequency Offset ' freq_fields{i} ' vs. k-mean \sigma_{l}'], direct, ...
        num_clusters, ['pkm_v_freq_' freq_fields{i}], 'log');
    
	relations.pkm.fft.(freq_fields{i}).R2 = R2;
 	relations.pkm.fft.(freq_fields{i}).R = R;
    relations.pkm.fft.(freq_fields{i}).p = p;
    
    % compare to gmm likelihood score
    [~, R2, R, p] = score_plot(full_name, sub_name, pgmm, fft.(freq_fields{i}), 'GMM \sigma_{l}', ...
        ['Frequency Offset ' freq_fields{i}], ['Frequency Offset ' freq_fields{i} ' vs. GMM \sigma_{l}'], direct, ...
        num_clusters, ['pgmm_v_freq_' freq_fields{i}], 'log');
    
	relations.pgmm.fft.(freq_fields{i}).R2 = R2;
 	relations.pgmm.fft.(freq_fields{i}).R = R;
    relations.pgmm.fft.(freq_fields{i}).p = p;
end
end

%% Phase Relations

if ~isempty(phase_comp)
for i = 1:length(phase_fields)
    
    % compare to k-means trasnition norm score
    [~, R2, R, p] = score_plot(full_name, sub_name, fkm, phase.(phase_fields{i}), 'k-means \sigma_{d}', ...
        ['Phase a_{1} a_{2} ' phase_fields{i}], ['Phase a_{1} a_{2} ' phase_fields{i} ' vs. k-mean \sigma_{d}'], direct, ...
        num_clusters, ['fkm_v_phase_' phase_fields{i}], 'linear');
    
	relations.fkm.phase.(phase_fields{i}).R2 = R2;
 	relations.fkm.phase.(phase_fields{i}).R = R;
    relations.fkm.phase.(phase_fields{i}).p = p;
    
    % compare to gmm trasnition norm score
    [~, R2, R, p] = score_plot(full_name, sub_name, fgmm, phase.(phase_fields{i}), 'GMM \sigma_{d}', ...
        ['Phase a_{1} a_{2} ' phase_fields{i}], ['Phase a_{1} a_{2} ' phase_fields{i} ' vs. GMM \sigma_{d}'], direct, ...
        num_clusters, ['fgmm_v_phase_' phase_fields{i}], 'linear');
    
	relations.fgmm.phase.(phase_fields{i}).R2 = R2;
 	relations.fgmm.phase.(phase_fields{i}).R = R;
    relations.fgmm.phase.(phase_fields{i}).p = p;
    
    % compare to k-mean likelihood score
    [~, R2, R, p] = score_plot(full_name, sub_name, pkm, phase.(phase_fields{i}), 'k-means \sigma_{l}', ...
        ['Phase a_{1} a_{2} ' phase_fields{i}], ['Phase a_{1} a_{2} ' phase_fields{i} ' vs. k-mean \sigma_{l}'], direct, ...
        num_clusters, ['pkm_v_phase_' phase_fields{i}], 'linear');
    
	relations.pkm.phase.(phase_fields{i}).R2 = R2;
 	relations.pkm.phase.(phase_fields{i}).R = R;
    relations.pkm.phase.(phase_fields{i}).p = p;
    
    % compare to gmm likelihood score
    [~, R2, R, p] = score_plot(full_name, sub_name, pgmm, phase.(phase_fields{i}), 'GMM \sigma_{l}', ...
        ['Phase a_{1} a_{2} ' phase_fields{i}], ['Phase a_{1} a_{2} ' phase_fields{i} ' vs. GMM \sigma_{l}'], direct, ...
        num_clusters, ['pgmm_v_phase_' phase_fields{i}], 'linear');
    
	relations.pgmm.phase.(phase_fields{i}).R2 = R2;
 	relations.pgmm.phase.(phase_fields{i}).R = R;
    relations.pgmm.phase.(phase_fields{i}).p = p;
end
end

end



    