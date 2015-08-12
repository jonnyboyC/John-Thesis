function frequency = compare_freq(frequency, modal_amp, sample_freq, completed, ...
    target_freq, model, sub_model, file)
% FREQ_SCORE attempt to match peaks in FRF that coorespond best to some
% target_frequency
%
%   frequency = FREQ_SCORE(frequency, modal_amp, sample_freq, target_freq,
%       model, sub_mode, file)


[~, ~, ~, loc] = calc_FRF2(modal_amp, sample_freq, 10);

for i = 1:3
    [best, distance] = knnsearch(loc(:,i), target_freq);
    frequency.(file.name).(model).(sub_model).(['amp' num2str(i) '_diff']) = distance;
    frequency.(file.name).(model).(sub_model).(['amp' num2str(i) '_val']) = loc(best,i);
    frequency.(file.name).(model).(sub_model).mode = i;
end
frequency.(file.name).(model).(sub_model).completed = completed.(model).(sub_model);
end