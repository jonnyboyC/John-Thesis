function TKE = calc_energy(TKE, integration, completed, modal_amp, model, sub_model, file, modes, MOD)
% CALC_ENERGY calculate system energy in score_all, packs results into TKE
% struct
%
%   TKE = CALC_ENERGY(TKE, integration, modal_amp, model, submodel, ...
%           file, modes)

% Pull off first mode if not modified
if ~MOD
    sim_modal_amp = integration.modal_amp.(model).(sub_model)(:,2:end);
else
    sim_modal_amp = integration.modal_amp.(model).(sub_model);
end

% Get appropraite modal_amps from data
temp_modal_amp = modal_amp(:,modes+1);

% Calculate or store mean and std of energy
TKE.(file.name).(model).(sub_model).mean_diff = abs(mean(sum(1/2*sim_modal_amp'.^2)) ...
    - mean(sum(1/2*temp_modal_amp'.^2)));
TKE.(file.name).(model).(sub_model).median_diff = abs(median(sum(1/2*sim_modal_amp'.^2)) ...
    - median(sum(1/2*temp_modal_amp'.^2)));
TKE.(file.name).(model).(sub_model).std_diff = abs(std(sum(1/2*sim_modal_amp'.^2)) ...
    - std(sum(1/2*temp_modal_amp'.^2)));
TKE.(file.name).(model).(sub_model).completed = completed.(model).(sub_model);

end