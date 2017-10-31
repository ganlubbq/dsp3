
clear

nsample = 200;
sim_tx = 1e-4;
wvar = 1e-2;
pnoise = phase_noise(nsample, sim_tx, 0);
wnoise = gaussian_noise(size(pnoise,1), size(pnoise,2), wvar, 'linear', 'real');
x = pnoise + wnoise;
% figure; plot(1:nsample, pnoise, 1:nsample, pnoise + wnoise); grid on
s_0 = 0;
M_0 = sim_tx;
for ii = 1 : nsample
    if ii == 1
        s(ii, 1) = s_0;
        M(ii, 1) = M_0 + sim_tx;
        K(ii) = M(ii, 1) / (M(ii, 1) + wvar);
        s(ii, 2) = s(ii, 1) + K(ii) * (x(ii) - s(ii, 1));
        M(ii, 2) = (1 - K(ii)) * M(ii, 1);
    else
        s(ii, 1) = s(ii - 1, 2);
        M(ii, 1) = M(ii - 1, 2) + sim_tx;
        K(ii) = M(ii, 1) / (M(ii, 1) + wvar);
        s(ii, 2) = s(ii, 1) + K(ii) * (x(ii) - s(ii, 1));
        M(ii, 2) = (1 - K(ii)) * M(ii, 1);
    end
end
figure; plot(1:nsample, x, 1:nsample, s(:,2)); grid on
figure; plot(1:nsample, pnoise, 1:nsample, s(:,2)); grid on



% mn = 4;
% a = constellation(mn);
% sp = sum(abs(a).^2) / mn;
% snr = 13;
% bit_tx = randi([0 1], 2, nsample);
% sim_tx = symbolizer_mqam(bit_tx);
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
% wnoise = gaussian_noise(size(sim_tx,1), size(sim_tx,2), sigma2, 'dbw', 'complex');
% sim_tx = sim_tx .* exp(1i * pnoise) + wnoise;
% pilot = pilot .* exp(1i * pnoise) + wnoise;
% 
% H = frequency_response(nsample, 1, 0.05, 0.01, 'rc');
% sim_rx = sim_tx .* conj(sign(ifft(fft(pilot).* H)));
% % scatterplot(sim_rx);
% figure; plot(unwrap(angle(sign(ifft(fft(pilot).* H))))); hold on; plot(pnoise);
% ber = nnz(bit_tx - slicer_mqam(sim_rx, mn)) / (nsample * 2)

