
clear

nsample = 2000;
% pvar = 1e-4;
% wvar = 1e-2;
% pnoise = phase_noise(nsample, pvar, 0);
% wnoise = gaussian_noise(size(pnoise,1), size(pnoise,2), wvar, 'linear', 'real');
% x = pnoise + wnoise;
% % figure; plot(1:nsample, pnoise, 1:nsample, pnoise + wnoise); grid on
% s_0 = 0;
% M_0 = pvar;
% for ii = 1 : nsample
%     if ii == 1
%         s(ii, 1) = s_0;
%         M(ii, 1) = M_0 + pvar;
%         K(ii) = M(ii, 1) / (M(ii, 1) + wvar);
%         s(ii, 2) = s(ii, 1) + K(ii) * (x(ii) - s(ii, 1));
%         M(ii, 2) = (1 - K(ii)) * M(ii, 1);
%     else
%         s(ii, 1) = s(ii - 1, 2);
%         M(ii, 1) = M(ii - 1, 2) + pvar;
%         K(ii) = M(ii, 1) / (M(ii, 1) + wvar);
%         s(ii, 2) = s(ii, 1) + K(ii) * (x(ii) - s(ii, 1));
%         M(ii, 2) = (1 - K(ii)) * M(ii, 1);
%     end
% end
% figure; plot(1:nsample, x, 1:nsample, s(:,2)); grid on
% figure; plot(1:nsample, pnoise, 1:nsample, s(:,2)); grid on



mn = 4;
a = constellation(mn);
sp = sum(abs(a).^2) / mn;
snr = 20;
bitTx = randi([0 1], 2, nsample);
symTx = symbolizer_mqam(bitTx);
symRef = symTx;

% define phase noies
txLaserPnVar = 1e-5;
phaseNoise = phase_noise(nsample, txLaserPnVar, 0);
phaseNoise = phaseNoise(:);

% debug, one could also test fixed phase error
% phaseNoise = 0;

% add noise
sigma2 = 10*log10(sp) - snr;
symTx = symTx .* exp(1i * phaseNoise) + gaussian_noise(size(symTx,1), size(symTx,2), sigma2, 'dbw', 'complex');

H = frequency_response(nsample, 1, 0.05, 0.001, 'rc');
symRec_1 = symTx .* exp(-1i * unwrap(angle(ifft(fft(symTx.^4).*H))/4)) .* exp(1i*pi/4);
symRec_2 = symTx .* exp(-1i * unwrap(angle(ifft(fft(symTx.^4).*H))/4)) .* exp(-1i*pi/4);
% scatterplot(symRec)

symRx = normalization(symRec_1, mn);
bitrx = slicer_mqam(symRx, mn);
ber_1 = nnz(bitTx - bitrx) / (nsample * 2);

symRx = normalization(symRec_2, mn);
bitrx = slicer_mqam(symRx, mn);
ber_2 = nnz(bitTx - bitrx) / (nsample * 2);

ber = [ber_1, ber_2]
