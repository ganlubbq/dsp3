
clear

nsample = 10000;
pvar = 1e-4;
wvar = 1e-2;
pnoise = phase_noise(nsample, pvar, 0);
wnoise = gaussian_noise(size(pnoise,1), size(pnoise,2), wvar, 'linear', 'real');
x = pnoise + wnoise;
% figure; plot(1:nsample, pnoise, 1:nsample, pnoise + wnoise); grid on
s_0 = 0;
M_0 = pvar;
for ii = 1 : nsample
    if ii == 1
        s(ii, 1) = s_0;
        M(ii, 1) = M_0 + pvar;
        K(ii) = M(ii, 1) / (M(ii, 1) + wvar);
        s(ii, 2) = s(ii, 1) + K(ii) * (x(ii) - s(ii, 1));
        M(ii, 2) = (1 - K(ii)) * M(ii, 1);
    else
        s(ii, 1) = s(ii - 1, 2);
        M(ii, 1) = M(ii - 1, 2) + pvar;
        K(ii) = M(ii, 1) / (M(ii, 1) + wvar);
        s(ii, 2) = s(ii, 1) + K(ii) * (x(ii) - s(ii, 1));
        M(ii, 2) = (1 - K(ii)) * M(ii, 1);
    end
end
% figure; plot(1:nsample, x, 1:nsample, s(:,2)); grid on
figure; plot(1:nsample, pnoise, 1:nsample, s(:,2)); grid on

x = exp(1i * pnoise) + wnoise;
fs = 2;
psd_x = spectrum_analyzer(x, fs, 'off');
psd_w = spectrum_analyzer(wnoise, fs, 'off');
freq = getFFTGrid(nsample, fs);
% pvar = 1e-5;
L = 4 * pvar ./ fs ./ (pvar^2 + 16 * pi * pi * freq.^2 ./ (fs)^2);
L = L / max(L) * max(psd_x);
% plot(fftshift(freq), dbw(fftshift(L)), 'LineWidth', 2);

% Wiener filter
H = L ./ (L + mean(psd_w));
y = ifft(fft(x) .* H);
figure; 
plot(unwrap(angle(x - wnoise))); hold on; 
plot(unwrap(angle(y))); grid on;


% mn = 4;
% a = constellation(mn);
% sp = sum(abs(a).^2) / mn;
% snr = 13;
% bit_tx = randi([0 1], 2, nsample);
% pvar = symbolizer_mqam(bit_tx);
% pilot = ones(nsample, 1);
% 
% % define phase noies
% pvar = 1e-5;
% pnoise = phase_noise(nsample, pvar, 0);
% pnoise = pnoise(:);
% 
% % debug, one could also test fixed phase error
% % pnoise = 0;
% 
% % add noise
% sigma2 = 10*log10(sp) - snr;
% wnoise = gaussian_noise(size(pvar,1), size(pvar,2), sigma2, 'dbw', 'complex');
% pvar = pvar .* exp(1i * pnoise) + wnoise;
% pilot = pilot .* exp(1i * pnoise) + wnoise;
% 
% H = frequency_response(nsample, 1, 0.05, 0.01, 'rc');
% sim_rx = pvar .* conj(sign(ifft(fft(pilot).* H)));
% % scatterplot(sim_rx);
% figure; plot(unwrap(angle(sign(ifft(fft(pilot).* H))))); hold on; plot(pnoise);
% ber = nnz(bit_tx - slicer_mqam(sim_rx, mn)) / (nsample * 2)

