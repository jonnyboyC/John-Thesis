function frequency = freq_score(frequency, modal_amp, sample_freq, target_freq, ...
    model, sub_model, file)
% FREQ_SCORE attempt to match peaks in FRF that coorespond best to some
% target_frequency
%
%   frequency = FREQ_SCORE(frequency, modal_amp, sample_freq, target_freq,
%       model, sub_mode, file)

frequency.(file.name).(model).(sub_model).diff = inf;
frequency.(file.name).(model).(sub_model).val = inf;
frequency.(file.name).(model).(sub_model).mode = 0;

[~, ~, peaks, loc] = calc_FRF2(modal_amp, sample_freq, 10);
for i = 1:size(peaks,1)
    [best, distance] = knnsearch(target_freq, loc);
    if distance < frequency.(file.name).(model).(sub_model).diff
        frequency.(file.name).(model).(sub_model).diff = distance;
        frequency.(file.name).(model).(sub_model).val = best;
        frequency.(file.name).(model).(sub_model).mode = i;
    end
end
end