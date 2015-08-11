function TKE = calc_energy(TKE, modal_amp_sim, completed, modal_amp, model, sub_model, file, modes)
% CALC_ENERGY calculate system energy in score_all, packs results into TKE
% struct
%
%   TKE = CALC_ENERGY(TKE, integration, modal_amp, model, submodel, ...
%           file, modes)

% Get appropraite modal_amps from data
temp_modal_amp = modal_amp(:,modes);

% Calculate or store mean and std of energy
TKE.(file.name).(model).(sub_model).mean_diff = abs(mean(sum(1/2*modal_amp_sim'.^2)) ...
    - mean(sum(1/2*temp_modal_amp'.^2)));
TKE.(file.name).(model).(sub_model).median_diff = abs(median(sum(1/2*modal_amp_sim'.^2)) ...
    - median(sum(1/2*temp_modal_amp'.^2)));
TKE.(file.name).(model).(sub_model).std_diff = abs(std(sum(1/2*modal_amp_sim'.^2)) ...
    - std(sum(1/2*temp_modal_amp'.^2)));
TKE.(file.name).(model).(sub_model).completed = completed.(model).(sub_model);

end