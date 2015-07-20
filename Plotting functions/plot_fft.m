function handle = plot_fft(modal_amp, num_modes, sample_freq, xlim, direct, type, custom, id)
% Calculate the modal frequency response using fft

if num_modes > 4
    num_plot = 1:4;
else
    num_plot = 1:num_modes;
end

% Number of points per windows, set to the next power of 2
window_samples = ceil(size(modal_amp,1)/4);
NFFT    = 2^nextpow2(window_samples); 
windows = floor(size(modal_amp, 1)/NFFT);   % Number of windows
T       = NFFT/sample_freq;                 % window sample by sampling rate
d_freq  = 1/T;                              % frequency resolution
fspan   = 0:d_freq:sample_freq-d_freq;

window  = hanning(NFFT);        % Hanning windows of width NFFT
start   = 1;                    % Start location  
finish  = NFFT;                 % End location 
freq_response = 0;              % frequency response

% Calculate fft for selected hanning window
for i = 1:windows
    modal_amp_win = modal_amp(start:finish,:).*(window*ones(1,num_modes));
    freq_response_temp  = fft(modal_amp_win);
    freq_response_temp = abs(freq_response_temp)/sample_freq;
    freq_response = freq_response + freq_response_temp;
    start = start + NFFT;
    finish = finish + NFFT;
end

freq_response_dB = 20*log10(freq_response);

% setup plot frequency response
hf = figure;
axf = newplot;

% Label axis and title
axf.XLabel.String = 'frequency (Hz)';
axf.YLabel.String = 'Frequency Response (dB)';
axf.Title.String  = ['Modal Frequency Response ' id];
axf.XLim = xlim;

% setup plot phase angle

% Return if we have no information to avoid throwing error
if size(freq_response_dB,1) >= 2
    % plot frequency response 
    plot(fspan, freq_response_dB(:, num_plot));
    
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