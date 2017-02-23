% According to parseval's theorem, the signal power in time domain equals
% the signal power in frequency domain. However, due to normalization used
% in Matlab, the power integrated in frequency domain via fft should be
% scaled by 1/N
%
%   power_in_time = sum(abs(fft(signal)).^2) / N = power_in_frequency
%
% Therefore, the power density displayed on spectrum analyzer should be
%
%   psd_on_sa = abs(fft(signal)).^2 / N
% 
% so that the integration of it would be the total power
%
%   total_power = sum(psd_on_sa)
%
% However, the psd is not power per Hz, rather power per frequency slot of
% spectrum analyzer, varying with the frequency resolution
%
% Copyright DEFAULT
%
function [PSD, freqVect] = spectrumAnalyzer(x, freqVect)

if nargin < 2
	nSamples = length(x);
	freqVect = [(0:nSamples/2-1)'; flipud(-(1:(nSamples/2))')] * 1.0 / nSamples;
else
    nSamples = length(x);
end

freqResolution = (max(freqVect) - min(freqVect)) / (nSamples - 1);

% power in one freq slot
PSD = abs(fft(x)) .^ 2 / nSamples;

figure(); grid on
plot(fftshift(freqVect), 10*log10(fftshift(PSD))); 

xlim([min(freqVect), max(freqVect)]);

% set sensitivity by limiting the y axis
% to do

xlabel('Frequency (Hz)'); 
ylabel('PSD (dB)');
legend(sprintf('Frequency Resolution %.2f MHz', freqResolution/1e6));

return

