function frequency = compare_freq(frequency, modal_amp, sample_freq, completed, ...
    target_freq, model, sub_model, file)
% FREQ_SCORE attempt to match peaks in FRF that coorespond best to some
% target_frequency
%
%   frequency = FREQ_SCORE(frequency, modal_amp, sample_freq, target_freq,
%       model, sub_mode, file)

frequency.(file.name).(model).(sub_model).diff = inf;
frequency.(file.name).(model).(sub_model).val = inf;
frequency.(file.name).(model).(sub_model).mode = 0;

[~, ~, peaks, loc] = calc_FRF2(modal_amp, sample_freq, 10);
for i = 1:size(peaks,2)
    [best, distance] = knnsearch(loc(:,i), target_freq);
    if distance < frequency.(file.name).(model).(sub_model).diff
        frequency.(file.name).(model).(sub_model).diff = distance;
        frequency.(file.name).(model).(sub_model).val = loc(best,i);
        frequency.(file.name).(model).(sub_model).mode = i;
    end
end
frequency.(file.name).(model).(sub_model).completed = completed.(model).(sub_model);
end