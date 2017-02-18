%% TEST SCRIPT FOR CALCULATING THEORETICAL BER VS OSNR
% CALCULATING THEORETICAL BER OF M-QAM SYSTEM WITH GAUSSIAN NOISE
% Denoted to pulse shaping, multiple samples per symbol, observing the ISI
% effet.
%%
function [] = TheoreticalBERv3(input_k)

if nargin<1;
    input_k = 1;
end

nSymbol = 2^16;

% NUMBER OF BIT PER SYMBOL
k = input_k;

refbit = randi([0 1],k,nSymbol);

% mapping bit to symbol
sym = symbolizerGrayQam(refbit);

% symbol power
ps = sum(abs(sym).^2)/nSymbol; 

snr = -10:0.5:10; % in dB

sps = 16;

nSamples = sps*nSymbol;

% Upsampling
if iscolumn(sym)
    sym = sym.';
end
sym = repmat(sym,sps,1);
sym = sym(:);

% get a freq domain raised cosine filter response
Rs = 1;
Fs = sps;
freqVect = getFFTGrid(nSamples,Fs);
alpha = 0.36; mode = 0;
H = calcRcosResponse(nSamples,Fs,Rs,alpha,mode);

% filtering signal in frequency domain
symi = real(ifft(fft(real(sym)).*H));
symq = real(ifft(fft(imag(sym)).*H));

sym = symi + 1j*symq;

ps = sum(abs(sym).^2)/nSamples; % symbol power

h = plotEyeDiagram(sym(1:8192),16,'e');



end

