function [psd, pxx, freq] = spectrum_analyzer(x, fs)
% According to parseval's theorem, the signal power in time domain equals
% the signal power in frequency domain. However, due to normalization used
% in Matlab, the power integrated in frequency domain via fft should be
% scaled by 1/N
%
%   power_in_time = sum(abs(signal).^2) / N = power_in_frequency =
%   sum(abs(fft(signal)).^2) / N^2
%
% Therefore, the power density displayed on spectrum analyzer should be
%
%   psd_on_sa = abs(fft(signal)).^2 / N^2
% 
% so that the integration of it would be the total power
%
%   total_power = sum(psd_on_sa)
%
% However, the psd is not power per Hz, rather power per frequency slot of
% spectrum analyzer, varying with the frequency resolution.
%
% To show the true PSD, the psd should be divided further by the frequency
% resolution.

if nargin < 2
    fs = 2.0;
end

nsample = length(x);
if mod(nsample, 2)
  freq = [(0 : 1 : (nsample - 1)/2), (-(nsample - 1)/2 : 1 : -1)] * fs / nsample;
else
  freq = [(0 : 1 : nsample/2 - 1), (-nsample/2 : 1 : -1)] * fs / nsample;
end

resolution = (max(freq) - min(freq)) / (nsample - 1);

% power in one freq slot
psd = abs(fft(x)) .^ 2 / (nsample * nsample);

% periodogram is the true power density
pxx = psd / (fs / nsample);

% figure(99);
plot(fftshift(freq), 10*log10(fftshift(psd)));
xlim([min(freq), max(freq)]);
ylim([10*log10(eps), 10]);
xlabel('Frequency (Hz)');
ylabel('Power per sample (dB)');
grid on;

if resolution >= 1e6
    title(sprintf('Spectrum Analyzer: Resolution %.2f MHz', resolution/1e6));
elseif resolution >= 1e3
    title(sprintf('Spectrum Analyzer: Resolution %.2f KHz', resolution/1e3));
else
    title(sprintf('Spectrum Analyzer: Resolution %.2f Hz', resolution));
end

return
