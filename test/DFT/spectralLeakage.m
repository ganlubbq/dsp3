% A block of data samples can be obtained by multiplying the infinite long
% input discrete signal by a rectangular window function, corresponding to
% convolute the Fourier transform of input siganl with a sinc function in
% frequency domain, which is the cause of so-called "spectral leakage".
% A Dirac delta in frequency domain will become the frequency response of
% the window function, hence the energy leaks into other frequencies.

clear
close all

%% global parameters
% number of DFT
N = 2048;
% freq of signal
f1 = 400;
f2 = 1000;
% sampling speed
fs = 8000;

% the maximal energy in frequency
fm1 = N * f1 / fs
fm2 = N * f2 / fs


%% signal with 2 distinct frequencies
signal_1 = cos(2 * pi * f1 * (0 : N-1) ./ fs);
signal_2 = cos(2 * pi * f2 * (0 : N-1) ./ fs);

% get the periodogram
temp = fft(signal_1 + signal_2);
psd = abs(temp(1 : N/2)).^2 / N /fs;

% however, due to finite resolution of DFT, the maximal energy of frequency
% that is not falling on the grid on DFT leaks into frequency bins around
% fm
% the spectral leakage is frequency dependent, for certain frequencies
% there is no spectral leakage as the fm coincides with one of the DFT
% bins
figure; plot(0 : N/2 - 1, dbw(psd)); grid on; legend('original DFT');


%% zero-padding the data samples will not increase the actual frequency
% resolution, but only to show more details of window pattern...even if
% frequency 2 falls exactly on one of the DFT grid, spectral leakage still
% can be observed
signal = [signal_1 + signal_2, zeros(1, 8000 - N)];
N = 8000;

% get the periodogram
temp = fft(signal);
psd = abs(temp(1 : N/2)).^2 / N / fs;

figure; plot(0 : N/2 - 1, dbw(psd)); grid on; legend('N = 2048, zero-padded to 8000');


%% increas nfft
N = 8000;

% the maximal energy in frequency
fm1 = N * f1 / fs
fm2 = N * f2 / fs

signal_1 = cos(2 * pi * f1 * (0 : N-1) ./ fs);
signal_2 = cos(2 * pi * f2 * (0 : N-1) ./ fs);

% get the periodogram
temp = fft(signal_1 + signal_2);
psd = abs(temp(1 : N/2)).^2 / N /fs;

figure; plot(0 : N/2 - 1, dbw(psd)); grid on; legend('NFFT = 8000');
