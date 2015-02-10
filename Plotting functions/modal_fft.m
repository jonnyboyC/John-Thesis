function handle = modal_fft(modal_amp, num_plot, window_samples, sample_freq, xlim, direct, MOD)
% Calculate the modal frequency response using fft

% Number of points per windows, set to the next power of 2
num_modes = size(modal_amp,2);
NFFT    = 2^nextpow2(window_samples); 
windows = floor(size(modal_amp, 1)/NFFT);   % Number of windows
T       = NFFT/sample_freq;                 % window sample by sampling rate
d_freq  = 1/T;                              % frequency resolution
fspan   = 0:d_freq:sample_freq-d_freq;

window  = hanning(NFFT);        % Hanning windows of width NFFT
start   = 1;                    % Start location  
finish  = NFFT;                 % End location 
freq_response = 0;              % frequency response
freq    = 0;                    % TODO

% Calculate fft for selected hanning window
for i = 1:windows
    modal_amp_win = modal_amp(start:finish,:).*(window*ones(1,num_modes));
    freq_response_temp  = fft(modal_amp_win);
    % TODO confirm we are dividing by num_elem
    freq_response_temp = abs(freq_response_temp)/sample_freq;
    freq_response = freq_response + freq_response_temp;
    start = start + NFFT;
    finish = finish + NFFT;
end

freq_response_dB = 20*log10(freq_response);

% Plot fft
h = figure;
ax = newplot;
if size(freq_response_dB,1) < 2
   if nargout == 1
       handle = h;
   end
   return
end
plot(fspan, freq_response_dB(:, num_plot));
ax.XLabel.String = 'frequency (Hz)';
ax.YLabel.String = 'Modal Amplitude';
ax.Title.String  = 'Modal Frequency Response';
ax.XLim = xlim;

leg_names = cell(size(num_plot,1),1);
for i = 1:size(num_plot,2)
   leg_names{i} = ['Modal FRF ' num2str(num_plot(i))]; 
end
legend(leg_names);

file_name = [direct '\Figures\Galerkin\FFT_'];
if nargin == 8
    file_name = [file_name MOD '_'];
end
file_name = [file_name num2str(size(modal_amp,2)) '_' num2str(ceil(sample_freq)) 'Hz'];
drawnow;

saveas(h, file_name, 'fig');

if nargout == 1
    handle = h;
end