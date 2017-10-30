mn = 4;
symlen = 10000;
a = constellation(mn);
sp = sum(abs(a).^2) / mn;
snr = 10;
bitTx = randi([0 1], bitpersym, symlen);
symTx = symbolizer_mqam(bitTx);
symRef = symTx;

% define phase noies
txLaserPnVar = 1e-4;
phaseNoise = phase_noise(symlen, txLaserPnVar, 0);
phaseNoise = phaseNoise(:);

% debug, one could also test fixed phase error
% phaseNoise = 0;

% add noise
sigma2 = 10*log10(sp) - snr;
symTx = symTx .* exp(1i * phaseNoise) + gaussian_noise(size(symTx,1), size(symTx,2), sigma2, 'dbw', 'complex');
