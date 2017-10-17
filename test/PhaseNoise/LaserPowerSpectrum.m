% Ref: 
%
% Magarini, Maurizio, et al. "Empirical modeling and simulation of
% phase noise in long-haul coherent optical transmission systems." Optics
% Express 19.23 (2011): 22455-22461.
%
% Di Domenico, Gianni, Stéphane Schilt, and Pierre Thomann. "Simple
% approach to the relation between laser frequency noise and laser line
% shape." Applied optics 49.25 (2010): 4801-4807.
%
% 
% Try sigma2 = 2. The spectral line shape looks gaussian.
% Try sigma2 = 2e-3. The spectral line shape looks lorentzian.
% Try even smaller sigma2 = 2e-6. Bad model.
%
% A noncausal Wiener filter is also demonstrated to estimate the
% target, i.e., pilot-tone with phase noise exp(i*theta), from its noisy
% measurement.
%
% Ref: MIT course 6.011 chapter 11 Wiener filter
%
% The Wiener is not practical as it requires precise knowledge of target
% PSD

clear
close all

fs = 2e6;
nsample = 10^5;
t = 0 : (1/fs) : (nsample-1)/fs;
freq = getFFTGrid(nsample, fs);

% sigma2 = 2;
sigma2 = 2e-4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1 : 2
    % random walk phase noise, the spectrum of random walk has Lorentzian shape
    pn = phase_noise(nsample, sigma2, 0);
    
    % WGN
    w = gaussian_noise(size(pn,1), size(pn,2), .1, 'linear', 'complex');
    
    % data model
    x = exp(1i * pn) + w;
    
    % the signal has another Lorentzian spectrum line shape
    psd(ii, :) = abs(fft(x)) .^ 2 / (nsample * nsample);
    psd_s(ii, :) = abs(fft(x - w)) .^ 2 / (nsample * nsample);
    psd_w(ii, :) = abs(fft(w)) .^ 2 / (nsample * nsample);
end

signalSpectrum = mean(abs(psd));
wgnSpectrum = mean(abs(psd_w));
laserPowerSpectrum = mean(abs(psd_s));

figure; plot(fftshift(freq), dbw(fftshift(wgnSpectrum))); grid on; box on; hold on
% plot the measured PSD
plot(fftshift(freq), dbw(fftshift(signalSpectrum)));
% plot the target PSD
plot(fftshift(freq), dbw(fftshift(laserPowerSpectrum)));

% theoretical model of target PSD, todo, need the coeff at the center
L = 4 * sigma2 ./ fs ./ (sigma2^2 + 16 * pi * pi * freq.^2 ./ (fs)^2);
% normalize
L = L / max(L) * max(laserPowerSpectrum);
plot(fftshift(freq), dbw(fftshift(L)), 'LineWidth', 2);

% Wiener filter
H = L ./ (L + mean(wgnSpectrum));
plot(fftshift(freq), dbw(fftshift(H)), 'LineWidth', 2);
legend('WGN', 'Measured PSD', 'Target PSD', 'Lorentzian', 'Wiener filter freq. resp.');

y = ifft(fft(x) .* H);
figure; plot(unwrap(angle(x - w))); grid on; hold on; plot(unwrap(angle(y)));
