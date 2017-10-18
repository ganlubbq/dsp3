

clear
close all


%% global parameters
% number of DFT
nfft = 4096;
% freq of signal
f1 = 1000;
% sampling speed
fs = 8000;
%
nsample = 3850;


%% signal with 2 distinct frequencies
signal_1 = cos(2 * pi * f1 * (0 : nsample - 1) ./ fs);

% get the periodogram, zero-padding the signal to nfft, power of 2
signal = [signal_1, zeros(1, nfft - nsample)];

%
freq = getFFTGrid(nfft, fs);
psd = spectrum_analyzer(signal, fs);


%%

