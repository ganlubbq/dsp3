% Main reference: Analog Devices Application Note AN-1067
clear 

fc = 1e6;
fs = 4e6;
nsample = 10^6;
t = 0 : (1/fs) : (nsample-1)/fs;

% sinusoid phase noise: the power of carrier is distributed to spikes at
% side-band, spacing as the frequency of phase noise and decaying as the
% increasing order of Bessel function of the first kind at the amplitude of
% phase noise

% the amplitude of the first spike aside the carrier is roughly half of the
% amplitude of sinusoid phase noise (dBc, with respect to the carrier), and
% all the other higher order spikes are negligible when the phase noise is
% small
pn = 0.1 * cos(2*pi*100e3*t);
% data model
x = cos(2*pi*fc*t + pn);
mf = max(abs(fft(x)).^2)/(nsample*nsample);
figure(90); hold on;
spectrumAnalyzer(x, fs);

% white gaussian phase noise
pn = genWGN(1, nsample, 0.01, 'linear', 'real');
% data model
x = cos(2*pi*fc*t + pn);
spectrumAnalyzer(x, fs);

% random walk phase noise
pn = genLaserPhaseNoise(nsample, 0.01, 0);
% data model
x = cos(2*pi*fc*t + pn);
spectrumAnalyzer(x, fs); 

hold off;
