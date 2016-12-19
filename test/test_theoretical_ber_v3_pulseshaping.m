%% TEST SCRIPT FOR CALCULATING THEORETICAL BER VS OSNR
% CALCULATING THEORETICAL BER OF M-QAM SYSTEM WITH GAUSSIAN NOISE
%%

clear all
% close all

nSymbol = 10^6;

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
sps = 8;
Fs = sps*rs;
nSamples = sps*nSymbol

% Upsampling
sym = repmat(sym,sps,[]);
sym = sym(:);

freqVector = getFFTGrid(nSamples,Fs);

alpha = 0.36;
H = calcRcosResponse(nSamples,Fs,rs,alpha,0);

symi = real(ifft(fft(real(sym)).*H));
symq = real(ifft(fft(imag(sym)).*H));

figure; plot(symi(1:1000)); grid on

ps = sum(abs(sym).^2)/nSymbol; % symbol power





%for ndx = 1:length(snr)
%    pn = ps/idbw(snr(ndx));
%    noise = genWGN(size(sym,1),size(sym,2),pn,'linear','complex');
%    signal = sym + noise;
%    if k==1
%        bit = slicerBPSK(signal);
%    else
%        bit = slicerGrayQam(signal,2^k);
%    end
%    ber(ndx) = nnz(bit(:)-refbit(:))/(k*nSymbol);
%    disp(sprintf('snr = %d, ber = %.2e',snr(ndx),ber(ndx)));
%end
%t_ber = T_BER_SNR_mQAM(snr,2^k);
%% t_ber = (1/k)*3/2*erfc(sqrt(k*0.1*(10.^(EbNo_db/10))));
%figure; plot(snr,log10(ber),'s-',snr,log10(t_ber),'k-'); grid on

