function [psd, freqVect] = spectrumAnalyzer(x, fs)
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
% spectrum analyzer, varying with the frequency resolution
nsample = length(x);
if nargin < 2
    fs = 1.0;
end
if mod(nsample, 2)
  freqVect = [(0 : 1 : (nsample-1)/2), (-(nsample-1)/2 : 1 : -1)] * fs / nsample;
else
  freqVect = [(0 : 1 : nsample/2-1), (-nsample/2 : 1 : -1)] * fs / nsample;
end

freqResolution = (max(freqVect) - min(freqVect)) / (nsample - 1);

% power in one freq slot
psd = abs(fft(x)) .^ 2 / (nsample * nsample);

% figure(99);
plot(fftshift(freqVect), 10*log10(fftshift(psd)));
xlim([min(freqVect), max(freqVect)]);
ylim([10*log10(eps), 10]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB)');
grid on;

if freqResolution > 1e6
    title(sprintf('Spectrum Analyzer: FR %.2f MHz', freqResolution/1e6));
elseif freqResolution > 1e3
    title(sprintf('Spectrum Analyzer: FR %.2f KHz', freqResolution/1e3));
else
    title(sprintf('Spectrum Analyzer: FR %.2f Hz', freqResolution));
end

return
