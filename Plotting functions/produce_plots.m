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
u_scale     = plot_data.u_scale;
l_scale     = plot_data.l_scale;
id          = plot_data.id;
type        = plot_data.type;
plot_type   = plot_data.plot_type;
x           = plot_data.x;
y           = plot_data.y;
Mod         = plot_data.Mod;
bnd_idx     = plot_data.bnd_idx;

t_scale = u_scale/l_scale;
sample_freq = sample_freq*t_scale;
t = t/t_scale;

% If using Mod add a mode zero to make it work like Galerkin
if Mod == true
    [modal_amp, ~, pod_ut, pod_vt] = add_mode_zero(modal_amp, 1, pod_ut, pod_vt, ...
            results.mean_u, results.mean_v);
end

% TODO significant overhaul to this function
if any(strcmp(plot_type, 'video'))
    plot_prediction(pod_ut, pod_vt, x, y, bnd_idx, modal_amp, t, dimensions, direct, id)
end

% Strip mean_u, mean_v
modal_amp = modal_amp(:,2:end);

% Plot modal amplitude of the response
if any(strcmp(plot_type, 'amp'))
    plot_amp(modal_amp, t, direct, type, id);
end


% Plot modal fft
if any(strcmp(plot_type, 'fft'))
    if num_modes > 4
        num2plot = 1:4;
    else
        num2plot = 1:num_modes;
    end
    if size(t,1) > 8192
        window_size = 4096;
    else
        window_size = ceil(size(t,1)/3);
    end
    modal_fft(modal_amp, num2plot, window_size, ...
        sample_freq, fft_window, direct, type, id);
end
