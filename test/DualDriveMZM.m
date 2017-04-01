%% TEST SCRIPT FOR CALCULATING THEORETICAL BER VS OSNR
% CALCULATING THEORETICAL BER OF M-QAM SYSTEM WITH GAUSSIAN NOISE

clear
% close all

nSymbol = 4096;

% NUMBER OF BIT PER SYMBOL
k = 2;

% SYMBOL RATE
Rs = 28e9;

% NOISE BANDWIDTH 
bn = 12.5e9;

refbit = randi([0 1],k,nSymbol);

% mapping bit to symbol
sym = symbolizerGrayQam(refbit);

sps = 64;
Fs = sps * Rs;
nSamples = sps * nSymbol;

% Upsampling
sym = repmat(sym,sps,[]);
sym = sym(:);

freqVector = getFFTGrid(nSamples,Fs);

alpha = 0.36;
H = calcRCFreqResponse(nSamples,Fs,Rs,alpha,0);

order = 4;
bandwidth = 32 * Rs;
H = calcBesselResponse(nSamples,Fs,order,bandwidth);

symi = real(ifft(fft(real(sym)) .* H));
symq = real(ifft(fft(imag(sym)) .* H));

%figure; plot(real(sym(1:256)));

phiI = symi*pi/2 + pi/2;
phiQ = (symq*pi/2 + pi/2) + pi/2;

txdual = exp(1i * phiI) + exp(1i * phiQ);

plot(txdual,'.'); grid on;

%print('dualarm_mzm_iqmod','-depsc');

