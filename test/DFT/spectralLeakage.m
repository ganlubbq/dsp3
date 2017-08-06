
clear
close all

%% global parameters
% number of DFT
N = 2048;
% freq of signal
f = 400;
% sampling speed
fs = 8000;

% the maximal energy in frequency
fm = N * f / fs

signal = exp(1i * 2 * pi * f * (0 : N-1) ./ fs);

temp = fft(signal);
psd = abs(temp(1 : N/2)).^2;

% however, due to finite resolution of DFT, the maximal energy in frequency
% leaks into frequency bins around fm
figure; plot(0 : N/2-1, log10(psd)); grid on

% the spectral leakage is frequency dependent, for certain frequencies
% there is no spectral leakage as the fm coincides with one of the DFT
% bins.