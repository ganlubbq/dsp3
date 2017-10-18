%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A block of data samples can be obtained by multiplying the infinite long
% input discrete signal by a rectangular window function, corresponding to
% convoluting the Fourier transform of input siganl with a sinc function,
% which causes a "blurred" DTFT and the so-called "spectral leakage". 
% 
% Due to the convolution, a Dirac delta in DTFT will be broadened into the
% frequency response of the window function and therefore the energy of
% single frequency component leaks into other frequencies. 
%
% DFT is a sampled version of the blurred version of DTFT. A longer time
% window will produce a less blurred DTFT and therefore increase the DFT
% sampling resolution. For a given time window, zero-padding the time
% sequence will produce more details of the blurred DTFT, but will not
% increase the effective frequency resolution. Therefore, zero-padding is
% equivalent to oversampling a given DTFT but can not make the DTFT less
% blurred.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all


%% global parameters
% number of DFT
N = 256;
% freq of signal
f1 = 40;
f2 = 100;
% sampling speed
fs = 800;

freq = getFFTGrid(N, fs);

% the maximal energy in frequency, fm1 is off grid, fm2 is on grid
fm1 = N * f1 / fs
fm2 = N * f2 / fs


%% signal off grid
% due to finite resolution of DFT, the maximal energy of frequency that is
% not falling on the grid on DFT leaks into frequency bins around fm
signal_1 = cos(2 * pi * f1 * (0 : N-1) ./ fs);

% get the periodogram
temp = fft(signal_1);
% psd = abs(temp(1 : N/2)).^2 / N /fs;
% to reveal more details of DTFT estimation only
psd = 2 * abs(temp).^2 / N / N;

figure(1); 
plot(fftshift(freq), dbw(fftshift(psd)), 'LineWidth', 2); 
hold on; grid on; 
xlabel('Frequency (Hz)'); ylabel('|DFT|^2 / N (dB)');


%% signal on grid
% however, the spectral leakage is frequency dependent, for certain
% frequencies there is no spectral leakage as the fm coincides with one of
% the DFT bins
signal_2 = cos(2 * pi * f2 * (0 : N-1) ./ fs);

% get the periodogram
temp = fft(signal_2);
% psd = abs(temp(1 : N/2)).^2 / N /fs;
% to reveal more details of DTFT estimation only
psd = 2 * abs(temp).^2 / N / N;

figure(2); 
plot(fftshift(freq), dbw(fftshift(psd)), 'LineWidth', 2); 
hold on; grid on; 
xlabel('Frequency (Hz)'); ylabel('|DFT|^2 / N (dB)');


%% zero-padding the data samples will NOT increase the actual frequency
% resolution, but only to show more details of window pattern...even if
% frequency 2 falls exactly on one of the DFT grid, spectral leakage still
% can be observed
L = 800;
freq = getFFTGrid(L, fs);
signal = [signal_1, zeros(1, L - N)];

% get the periodogram
temp = fft(signal);
% psd = abs(temp(1 : L/2)).^2 / N /fs;
% to reveal more details of DTFT estimation only
psd = 2 * abs(temp).^2 / N / N;

figure(1); 
plot(fftshift(freq), dbw(fftshift(psd)), 'LineWidth', 1); 
xlabel('Frequency (Hz)'); ylabel('|DFT|^2 / N (dB)');
legend('nfft = 256', 'nfft = 256, zero-padding = 800');
xlim([0, max(freq)]);

L = 800;
signal = [signal_2, zeros(1, L - N)];

% get the periodogram
temp = fft(signal);
% psd = abs(temp(1 : L/2)).^2 / N /fs;
% to reveal more details of DTFT estimation only
psd = 2 * abs(temp).^2 / N / N;

figure(2); 
plot(fftshift(freq), dbw(fftshift(psd)), 'LineWidth', 1);
xlabel('Frequency (Hz)'); ylabel('|DFT|^2 / N (dB)');
legend('nfft = 256', 'nfft = 256, zero-padding = 800');
xlim([0, max(freq)]);
ylim([-60, 0]);


%% increase nfft will increase the actual frequency resolution
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
