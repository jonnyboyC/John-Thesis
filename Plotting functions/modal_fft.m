function handle = modal_fft(modal_amp, num_plot, sample_freq)
samples = size(modal_amp,1);
NFFT = 2^nextpow2(samples);
freq = sample_freq/2*linspace(0,1,NFFT/2+1);
freq_response = fft(modal_amp, NFFT, 1);
freq_response = freq_response(1:ceil((NFFT+1)/2), :)/NFFT;

ax = newplot;
h = semilogy(ax, freq', abs(freq_response(:, num_plot)));
ax.XLabel.String = 'frequency (Hz)';
ax.YLabel.String = 'Modal Amplitude';
ax.Title.String  = 'Modal Frequency Response';

leg_names = cell(size(num_plot,1),1);
for i = 1:size(num_plot,2)
   leg_names{i} = ['Modal FRF ' num2str(num_plot(i))]; 
end
legend(leg_names);


if nargout == 1
    handle = h;
end