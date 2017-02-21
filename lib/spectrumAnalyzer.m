function [PSD, freqVect] = spectrumAnalyzer(x, freqVect)
%
% Copyright DEFAULT
%

if nargin < 2
	nSamples = length(x);
	freqVect = [(0:nSamples/2-1)'; flipud(-(1:(nSamples/2))')] * 1.0 / nSamples;
else
    nSamples = length(x);
end

PSD = abs(fft(x) ./ nSamples) .^ 2;

figure(); grid on
plot(fftshift(freqVect), 10*log10(fftshift(PSD))); 

xlim([min(freqVect), max(freqVect)]);

% set sensitivity by limiting the y axis
% to do

xlabel('Frequency (Hz)'); 
ylabel('PSD (dB/Hz)');

return

