%% TEST SCRIPT FOR CALCULATING THEORETICAL BER VS OSNR
% CALCULATING THEORETICAL BER OF M-QAM SYSTEM WITH GAUSSIAN NOISE
%%
clear all
% close all

nSymbol = 1024;

% NUMBER OF BIT PER SYMBOL
k = 2

% SYMBOL RATE
rs = 28e9;

% NOISE BANDWIDTH 
bn = 12.5e9;

refbit = randi([0 1],k,nSymbol);

% mapping bit to symbol
if k==1
    sym = symbolizerBPSK(refbit);
else
    sym = symbolizerGrayQam(refbit);
end

snr = 0:0.5:10; % in dB
cr = 1; % coding rate
sps = 64;
Fs = sps*rs;
nSamples = sps*nSymbol

% Upsampling
sym = repmat(sym,sps,[]);
sym = sym(:);

freqVector = getFFTGrid(nSamples,Fs);

alpha = 0.36;
H = calcRCFreqResponse(nSamples,Fs,rs,alpha,0);
order = 4;
bandwidth = 0.75*rs;
H = calcBesselResponse(nSamples,Fs,order,bandwidth);

symi = real(ifft(fft(real(sym)).*H));
symq = real(ifft(fft(imag(sym)).*H));

%figure; plot(real(sym(1:256)));

phiI = symi*pi/2 + pi/2;
phiQ = (symq*pi/2 + pi/2) + pi/2;

txdual = exp(1j*phiI) + exp(1j*phiQ);

plot(txdual,'.'); grid on;

%print('dualarm_mzm_iqmod','-depsc');



