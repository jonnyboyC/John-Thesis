function [freq_response, fspan] = calc_FRF(modal_amp, num_modes, sample_freq)


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

end