%% Noise analysis on phase recovery based on pilot-tone
% Phase estimation algorithm based on pilot-tone contains 2 steps: LPF and
% phase rotation. Everything in the LPF range will therefore reduce to
% one-dimension and be positive only. It is easy to verify that the total
% power doesn't change but more power is shifted to the zero-frequency and
% the PSD level of WGN after phase recovery actually drops.

clear
close all

fc = 1e6;
fs = 20e6;
nsample = 10^7;
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
% domain, we can observe that everything including the noise is REDUCED TO
% ONE DIMENSION, while the seperated WGN remains two dimension after phase
% recovery.
figure(1); spectrumAnalyzer(xc, fs);
legend('Seperated WGN', 'Seperated WGN with lower PSD', 'Pilot with WGN with lower PSD');
% filter out the signal in the filtering zone
xcf = ifft(fft(xc) .* H.');
scatterplot(xcf); grid on; title('Pilot with WGN after phase recovery');
scatterplot(nc); grid on; title('Seperated WGN after phase recovery');

% interestingly, the (phase) noise power remains in the pilot is nonzero
% after phase recovery, the power of seperated WGN remains the same,
% however the combined noise power in the filtering zone is reduced after
% phase recovery, indicating that although the phase recovery is a linear
% process to the data mode (d * exp(ip) + w), the residual noise in the data
% and WGN part is correlated.
pc = exp(1i*pn) .* conj(xf) ./ abs(xf);
scatterplot(pc); grid on; title('Seperated pilot after phase recovery');

% so, if thers is another signal aside the pilot-tone waiting for the phase
% noise compensation, it's better not to be too close to the pilot-tone
% otherwise part of the singal will be REDUCED TO ONE DIMENSION as well.
% Although one can observe a noise power reduction in the signal......

