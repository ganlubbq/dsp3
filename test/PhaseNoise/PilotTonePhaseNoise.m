%% Try out different lowpass filtering applied to the pilot-tone with phase
% noise. Observe the spectrum line shape and the remaining phase noise
clear
% close all

fc = 1e6;
fs = 20e6;
nsample = 10^6;
t = 0 : (1/fs) : (nsample-1)/fs;

%% random walk phase noise
pn = genLaserPhaseNoise(nsample, 2*pi*1E-3, 0);
% a white gaussian additive noise
an = genWGN(1, nsample, .3, 'linear', 'complex');
% data model with zero freq pilot tone
x = exp(1i * pn) + an;

%% noise analysis
% in this section, i want to show that the phase recovery process based on
% a pilot-tone will lower the PSD of an seperated WGN in the filtering zone
% by shifting energy to the zero frequency component

% raised cosine
H = calcFilterFreqResp(nsample, fs, 0.01, 10e6, 'rc');
xf = ifft(fft(x) .* H.');
xc = x .* conj(xf) ./ abs(xf);
nc = an .* conj(xf) ./ abs(xf);

% observe there is a peak in the noise after compensation
figure(1); spectrumAnalyzer(an, fs); hold on
figure(1); spectrumAnalyzer(nc, fs);

% in this section, i want to show that the phase recovery process based on
% a pilot-tone will lower even more of the PSD of WGN added to that
% pilot-tone in the filtering zone by shifting even more energy to the zero
% frequency component

% now the noise peak at zero frequency is hidden in the peak of pilot-tone
% and therefore is not observable. However, by plotting the signal in time
% domain, we can observe that everything including the noise is reduced to
% one dimension, while the seperated WGN remains two dimension after phase
% recovery.
figure(1); spectrumAnalyzer(xc, fs);
legend('seperated WGN', 'seperated WGN with lower PSD', 'pilot with WGN with lower PSD');
% filter out the signal in the filtering zone
xcf = ifft(fft(xc) .* H.');
scatterplot(nc); grid on; title('seperated WGN after phase recovery');
scatterplot(xcf); grid on; title('pilot with WGN after phase recovery');

%% filter the pilot tone
% moving average
% ntaps = 5;
% taps = ones(1, ntaps) / ntaps;
% xf = filter(taps, 1, x);
% pf = filter(taps, 1, exp(1i*pn));
% anf = filter(taps, 1, an);

% raised cosine
H = calcFilterFreqResp(nsample, fs, 0.01, 20e6, 'rc');
xf = ifft(fft(x) .* H.');
pf = ifft(fft(exp(1i*pn)) .* H.');
anf = ifft(fft(an) .* H.');

%% gaussian snr with signal power of 1
snr = dbw(1 / calcrms(anf).^2);
fprintf('SNR is %.2f dB \n', snr);

%% remove the phase noise
xc = x .* conj(xf) ./ abs(xf);
pn_est = calcrms(xc - mean(xc)).^2;
fprintf('estimated power of noise is %.4g \n', pn_est)

pc = exp(1i*pn) .* conj(xf) ./ abs(xf);
nc = an .* conj(xf) ./ abs(xf);

figure(1); spectrumAnalyzer(nc, fs); hold on
figure(1); spectrumAnalyzer(xc, fs);
figure(1); spectrumAnalyzer(pc, fs); hold off

% figure(99); clf; hold on
% psd = spectrumAnalyzer(x, fs);
% psd = spectrumAnalyzer(xf, fs); 
% box on; hold off 

% scatterplot(x)
% scatterplot(xc)

figure(10); clf; hold on
plot(pn(1:150), 'LineWidth', 2);
plot(unwrap(angle(xf(1:150)))); 
legend('Actual phase noise', 'Estimated phase noise');
grid on; box on;

