function phase_comp = phase_comparision(phase_comp, phase, modal_amp_sim, ...
                    modal_amp, modes, completed, km_steady, gmm_steady, model, sub_model, file)
% PHASE_COMPARISION see if various phase relations hold in simulated model
%
%   phase_comp = PHASE_COMPARISION(phase_comp, phase, modal_amp_sim,
%       modal_amp, model, sub_model, files, modes)

if phase(1) > length(modes) || phase(2) > length(modes)
   fprintf('Phase modes not included in mode');
   return;
end

% emp_fft1 = fft(hanning(size(modal_amp,1)).*modal_amp(:,modes(phase(1))));
% emp_fft2 = fft(hanning(size(modal_amp,1)).*modal_amp(:,modes(phase(2))));

sim_fft1 = fft(hanning(size(modal_amp_sim,1)).*modal_amp_sim(:,phase(1)));
sim_fft2 = fft(hanning(size(modal_amp_sim,1)).*modal_amp_sim(:,phase(2)));

% emp_ang1 = angle(emp_fft1(2:end));
% emp_ang2 = angle(emp_fft2(2:end));

sim_ang1 = angle(sim_fft1(2:end));
sim_ang2 = angle(sim_fft2(2:end));

% emp_diff = abs(emp_ang1 - emp_ang2);
sim_diff = abs(sim_ang1 - sim_ang2);

phase_comp.(file.name).(model).(sub_model).mean_diff = ...
    abs(pi - mean(sim_diff));
phase_comp.(file.name).(model).(sub_model).std_diff = ...
    abs(std(sim_diff));
phase_comp.(file.name).(model).(sub_model).median_diff = ...
    abs(pi - median(sim_diff));
phase_comp.(file.name).(model).(sub_model).completed = ...
    completed.(model).(sub_model);
phase_comp.(file.name).(model).(sub_model).steady = ...
    km_steady.(model).(sub_model);
   
end