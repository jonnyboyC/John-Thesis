function cluster_plots(clusters_info)
% Will be garbage code :(

results_scores  = clusters_info.results_scores;
TKE             = clusters_info.TKE;
frequency       = clusters_info.frequency;
amplitude1      = clusters_info.amplitude1;
amplitude2      = clusters_info.amplitude2;
num_clusters    = clusters_info.num_clusters;
direct          = clusters_info.direct;

[fkm, full_name]   = merge_struct(results_scores, 'frob_km');
fgmm            = merge_struct(results_scores, 'frob_gmm');
pkm             = merge_struct(results_scores, 'like_km');
pgmm             = merge_struct(results_scores, 'like_gmm');

if isempty(full_name)
    return;
end

tke_mean        = merge_struct(TKE, 'mean_diff');
tke_std         = merge_struct(TKE, 'std_diff');
tke_median      = merge_struct(TKE, 'median_diff');

if ~isempty(frequency)
    fft_diff        = merge_struct(frequency, 'diff');
    fft_mode        = merge_struct(frequency, 'mode');
end

amp1_mean        = merge_struct(amplitude1, 'mean_diff');
amp1_std         = merge_struct(amplitude1, 'std_diff');
amp1_median      = merge_struct(amplitude1, 'median_diff');

amp2_mean        = merge_struct(amplitude2, 'mean_diff');
amp2_std         = merge_struct(amplitude2, 'std_diff');
amp2_median      = merge_struct(amplitude2, 'median_diff');

sub_name = cell(size(full_name));
for i = 1:length(full_name)
    underscore = strfind(full_name{i}, '_');
    sub_name{i} = full_name{i}(underscore+1:end);
end



h = score_plot(full_name, sub_name, pgmm, tke_median, 'GMM \sigma_{l}', 'TKE median', ...
    'test', direct, num_clusters, 'tester_tester');
end



    