function produce_plots(plot_data)
% Produce plots for Galerkin_Proj and MOD_Pod based on requests passed in
% plot_type
modal_amp   = plot_data.modal_amp; 
num_modes   = plot_data.num_modes;
t           = plot_data.t;
direct      = plot_data.direct;
pod_Ut      = plot_data.pod_Ut;
dimensions  = plot_data.dimensions;
fft_window  = plot_data.fft_window;
sample_freq = plot_data.sample_freq;
id          = plot_data.id;
type        = plot_data.type;
plot_type   = plot_data.plot_type;
X           = plot_data.X;
Mod         = plot_data.Mod;
bnd_idx     = plot_data.bnd_idx;
custom      = plot_data.custom;

% If using Mod add a mode zero to make it work like Galerkin
if Mod == true
    [modal_amp, pod_Ut] = add_mode_zero_mod(modal_amp, pod_Ut, plot_data.mean_U);
end

% decide if streamlines should be used
if any(strcmp(plot_type, 'video stream'))
    use_stream = true;
else
    use_stream = false;
end

% Create video reconstruction of the flow
if any(strcmp(plot_type, 'video') | strcmp(plot_type, 'video stream'))
    plot_prediction(pod_Ut, X, bnd_idx, modal_amp, t, dimensions, use_stream, direct, custom, id)
end

% Strip mean mode off
modal_amp = modal_amp(:,2:end);

% Plot modal amplitude of the response
if any(strcmp(plot_type, 'amp'))
    plot_amp(modal_amp, t, direct, type, custom, id);
end

% Plot system energy
if any(strcmp(plot_type, 'energy'))
    plot_energy(modal_amp, t, id, direct, type, custom);
end

% Plot system fft
if any(strcmp(plot_type, 'fft'))
    plot_fft(modal_amp, num_modes,...
        sample_freq, fft_window, direct, type, custom, id);
end
