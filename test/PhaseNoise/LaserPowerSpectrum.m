% Ref: 
% Magarini, Maurizio, et al. "Empirical modeling and simulation of
% phase noise in long-haul coherent optical transmission systems." Optics
% Express 19.23 (2011): 22455-22461.

clear
close all

fc = 1e6;
fs = 4e6;
nsample = 10^5;
t = 0 : (1/fs) : (nsample-1)/fs;
freq = getFFTGrid(nsample, fs);

sigma2 = 1e-2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1 : 100
    % random walk phase noise, the spectrum of random walk has Lorentzian shape
    pn = genLaserPhaseNoise(nsample, sigma2, 0);
    % data model
    x = exp(1i * pn);
    % the signal has another Lorentzian spectrum line shape
    psd(ii, :) = abs(fft(x)) .^ 2 / (nsample * nsample);
end
laserPowerSpectrum = mean(abs(psd));
figure; plot(fftshift(freq), dbw(fftshift(laserPowerSpectrum))); grid on; box on; hold on

% theoretical model, todo, need the coeff at the center
L = 4 * sigma2 ./ fs ./ (sigma2^2 + 16 * pi * pi * freq.^2 ./ (fs)^2);
plot(fftshift(freq), dbw(fftshift(L/max(L)*max(laserPowerSpectrum))), 'LineWidth', 2);


