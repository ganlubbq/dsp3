% Main reference: Analog Devices Application Note AN-1067
clear
close all

fc = 1e6;
fs = 4e6;
nsample = 10^5;
t = 0 : (1/fs) : (nsample-1)/fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sinusoid phase noise: the power of carrier is distributed to spikes at
% side-band, spacing as the frequency of phase noise and decaying as the
% increasing order of Bessel function of the first kind at the amplitude of
% phase noise

% the amplitude of the first spike aside the carrier is roughly half of the
% amplitude of sinusoid phase noise (dBc, with respect to the carrier), and
% all the other higher order spikes are negligible when the phase noise is
% small. In PSD, the 1st spike should be 2 * dbw(an / 2) lower than the
% carrier.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
an = 0.1;
pn = an * cos(2 * pi * 100e3 * t);
% data model
x = 2 * cos(2 * pi * fc * t + pn);
mf = max(abs(fft(x)).^2) / (nsample * nsample);
figure(90); hold on;
spectrumAnalyzer(x, fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% white gaussian phase noise
pn = genWGN(1, nsample, 0.02, 'linear', 'real');
% data model
x = 2 * cos(2 * pi * fc * t + pn);
% the signal has white spectrum outside the carrier
psd = spectrumAnalyzer(x, fs);
% sum all the freq. components from 0 to fs/2, excluding the carrier, will
% give the power of phase noise
prms = sum(psd(1:nsample/2)) - psd(nsample/4 + 1);
fprintf('estimated power of phase noise %.4f\n', prms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% random walk phase noise, the spectrum of random walk has Lorentzian shape
pn = genLaserPhaseNoise(nsample, 1e-4, 0);
% data model
x = 2 * cos(2 * pi * fc * t + pn);
% the signal has another Lorentzian spectrum line shape
psd = spectrumAnalyzer(x, fs); 
% summing up does not work anymore!
prms = sum(psd(1:nsample/2)) - psd(nsample/4 + 1);
fprintf('estimated power of phase noise %.4f\n', prms);
box on; hold off;
