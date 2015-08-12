function [amplitude1, amplitude2] = compare_amp(amplitude1, amplitude2, modal_amp_sim, ...
    modal_amp, modes, completed, km_steady, gmm_steady, model, sub_model, file)
% COMPARE_AMP compare predicted mean and std deviation of individual modal
% amplitudes
% 
%   amplitude = COMPARE_AMP(amplitude, modal_amp_sim, modal_amp, model,
%       sub_model, file, modes, completed)

amplitude1.(file.name).(model).(sub_model).mean_diff = ...
    abs(mean(modal_amp_sim(:,1)) - mean(modal_amp(:,modes(1))));
amplitude1.(file.name).(model).(sub_model).median_diff = ...
    abs(median(modal_amp_sim(:,1)) - median(modal_amp(:,modes(1))));
amplitude1.(file.name).(model).(sub_model).std_diff = ...
    abs(std(modal_amp_sim(:,1)) - std(modal_amp(:,modes(1))));
amplitude1.(file.name).(model).(sub_model).completed = completed.(model).(sub_model);
amplitude1.(file.name).(model).(sub_model).steady = km_steady.(model).(sub_model);


amplitude2.(file.name).(model).(sub_model).mean_diff = ...
    abs(mean(modal_amp_sim(:,2)) - mean(modal_amp(:,modes(2))));
amplitude2.(file.name).(model).(sub_model).median_diff = ...
    abs(median(modal_amp_sim(:,2)) - median(modal_amp(:,modes(2))));
amplitude2.(file.name).(model).(sub_model).std_diff = ...
    abs(std(modal_amp_sim(:,2)) - std(modal_amp(:,modes(2))));
amplitude2.(file.name).(model).(sub_model).completed = completed.(model).(sub_model);
amplitude2.(file.name).(model).(sub_model).steady = km_steady.(model).(sub_model);



