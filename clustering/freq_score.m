function frequency = freq_score(frequency, modal_amp, sample_freq, target_freq, ...
    model, sub_model, file)

[~, fspan, peaks, loc] = calc_FRF2(modal_amp, sample_freq, 10);

end