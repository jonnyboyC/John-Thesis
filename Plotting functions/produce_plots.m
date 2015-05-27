function produce_plots(plot_data)
% Produce plots for Galerkin_Proj and MOD_Pod based on requests passed in
% plot_type
modal_amp   = plot_data.modal_amp; 
num_modes   = plot_data.num_modes;
t           = plot_data.t;
direct      = plot_data.direct;
pod_ut      = plot_data.pod_ut;
pod_vt      = plot_data.pod_vt;
dimensions  = plot_data.dimensions;
fft_window  = plot_data.fft_window;
sample_freq = plot_data.sample_freq;
id          = plot_data.id;
type        = plot_data.type;
plot_type   = plot_data.plot_type;
x           = plot_data.x;
y           = plot_data.y;
Mod         = plot_data.Mod;
bnd_idx     = plot_data.bnd_idx;
custom      = plot_data.custom;

use_stream = false;

% If using Mod add a mode zero to make it work like Galerkin
if Mod == true
    [modal_amp, pod_ut, pod_vt] = add_mode_zero_mod(modal_amp, pod_ut, pod_vt, ...
            plot_data.mean_u, plot_data.mean_v);
end

% TODO significant overhaul to this function
if any(strcmp(plot_type, 'video stream'))
    use_stream = true;
end

if any(strcmp(plot_type, 'video') | strcmp(plot_type, 'video stream'))
    plot_prediction(pod_ut, pod_vt, x, y, bnd_idx, modal_amp, t, dimensions, use_stream, direct, custom, id)
end

% Strip mean_u, mean_v
modal_amp = modal_amp(:,2:end);

% Plot modal amplitude of the response
if any(strcmp(plot_type, 'amp'))
    plot_amp(modal_amp, t, direct, type, custom, id);
end

if any(strcmp(plot_type, 'energy'))
    plot_energy(modal_amp, t, id, direct, type, custom);
end


% Plot modal fft
if any(strcmp(plot_type, 'fft'))
    if num_modes > 4
        num2plot = 1:4;
    else
        num2plot = 1:num_modes;
    end
    if size(t,1) > 8192
        window_size = ceil(size(t,1)/4);
%         window_size = 8192;
    else
        window_size = ceil(size(t,1)/4);
    end
    modal_fft(modal_amp, num2plot, window_size, ...
        sample_freq, fft_window, direct, type, custom, id);
end
