function [PSD freqVect] = spectrumAnalyzer(x, freqVect)
%
% Copyright default
%

if nargin<2
	nSamples = length(x);
	freqVect = [(0:nSamples/2-1)'; flipud(-(1:(nSamples/2))')] *1.0/nSamples;
end

PSD = abs(fft(x)./nSamples).^2;

figure();
plot(fftshift(freqVect), 10*log10(fftshift(PSD))); grid on
xlim([min(freqVect) max(freqVect)]);
% set sensitivity by limiting the y axis
% to do

xlabel('Frequency (Hz)'); 
ylabel('PSD (dB/Hz)');

return

