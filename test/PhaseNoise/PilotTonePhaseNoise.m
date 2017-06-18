% Pilot-tone with random walk phase noise
clear
% close all

fc = 1e6;
fs = 4e6;
nsample = 10^6;
t = 0 : (1/fs) : (nsample-1)/fs;

% random walk phase noise
pn = genLaserPhaseNoise(nsample, 2*pi*10e3/10e6 * 0, 0);
% a white gaussian additive noise
an = genWGN(1, nsample, 1e-1, 'linear', 'complex');
% data model with zero freq pilot tone
x = exp(1i * pn) + an;

% filter the pilot tone
ntaps = 1;
taps = ones(1, ntaps) / ntaps;
xf = filter(taps, 1, x);
anf = filter(taps, 1, an);

% gaussian snr with signal power of 1
snr = dbw(1 / calcrms(anf).^2);
fprintf('SNR is %.2f dB \n', snr);

% remove the phase noise
xc = x .* conj(xf);
pn_est = calcrms(xc - mean(xc)).^2;
fprintf('estimated power of noise is %.4g \n', pn_est)

figure(99); clf; hold on
psd = spectrumAnalyzer(x, fs);
psd = spectrumAnalyzer(xf, fs); hold off 

% scatterplot(x)
% scatterplot(xc)

figure(10); clf; hold on
plot(pn(1:150));
plot(unwrap(angle(xf(1:150)))); grid on; hold off

