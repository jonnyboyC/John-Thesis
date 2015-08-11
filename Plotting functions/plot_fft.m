function handle = plot_fft(modal_amp, num_modes, sample_freq, xlim, direct, type, custom, id, varargin)
% Calculate the modal frequency response using fft

if num_modes > 8
    num_plot = 1:8;
else
    num_plot = 1:num_modes;
end

if nargin == 9
    plot_peaks = varargin{1};
end

% setup plot handles
hf = figure;
axf = newplot;

% Calculate frequency response
if plot_peaks
    [freq_response, fspan, peaks, loc] = calc_FRF2(modal_amp, sample_freq);
else
    [freq_response, fspan] = calc_FRF2(modal_amp, sample_freq);
end
% [freq_response, fspan] = calc_FRF(modal_amp, num_modes, sample_freq);

freq_response_dB = 20*log10(freq_response);

% setup plot phase angle

% Return if we have no information to avoid throwing error
if size(freq_response_dB,1) >= 2
    % plot frequency response 
    plot(axf, fspan, freq_response_dB(:, num_plot));
    
    if plot_peaks
        hold(axf, 'on');
        plot(axf, fspan(loc(:, num_plot)), 20*log10(peaks(:, num_plot)));
    end
    
    % Label axis and title
    axf.XLabel.String = 'frequency (Hz)';
    axf.YLabel.String = 'Frequency Response (dB)';
    axf.Title.String  = ['Modal Frequency Response ' id];
    axf.XLim = xlim;
    
    % Generate legend
    leg_names = cell(size(num_plot,1),1);
    for i = 1:size(num_plot,2)
        leg_names{i} = ['Modal FRF ' num2str(num_plot(i))];
    end
    legend(leg_names);
end


if custom
    direct_ext = [direct filesep 'Figures' filesep type filesep 'modes_' ...
        num2str(num_modes) '_custom' filesep 'fourier'];
else
    direct_ext = [direct filesep 'Figures' filesep type filesep 'modes_' ...
        num2str(num_modes) filesep 'fourier'];
end

if ~exist(direct_ext, 'dir') 
    mkdir(direct_ext);
end

file_name = [direct_ext filesep 'FFT_' strrep(id, ' ', '_')];
file_name = [file_name '_' num2str(ceil(sample_freq)) 'Hz'];
drawnow;

saveas(hf, file_name, 'fig');
saveas(hf, file_name, 'png')

if nargout == 1
    handle = hf;
end