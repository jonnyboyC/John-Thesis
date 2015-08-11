function [FRF, fspan, peaks, loc] = calc_FRF2(modal_amp, sample_freq, varargin)
samples = size(modal_amp,1);
modes = size(modal_amp,2);
NFFT = 2^nextpow2(samples);
window = hanning(samples);
[FRF, fspan] = periodogram(modal_amp, window, NFFT, sample_freq, 'power');

if nargin == 3
    pks = varargin{1};
else
    pks = 20;
end

if nargout >= 3 
    count = find(fspan > 5, 1,'first');
    peaks = zeros(pks, modes);
    loc = zeros(pks,modes);
    for i = 1:modes
        if length(FRF(:,1)) < pks
            break;
        end
        [peaks_temp, location] = findpeaks(FRF(:,i), 'SortStr', 'descend', ...
            'MinPeakDistance', count, 'NPeaks', pks);
        found_peaks = length(peaks_temp);
        loc(1:found_peaks,i) = fspan(location);
        peaks(1:found_peaks,i) = peaks_temp;
    end
end