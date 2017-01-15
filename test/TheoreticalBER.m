%% TEST SCRIPT FOR CALCULATING THEORETICAL BER VS OSNR
% CALCULATING THEORETICAL BER OF M-QAM SYSTEM WITH GAUSSIAN NOISE
%% QPSK
clear

nSymbol = 10^6;

% NUMBER OF BIT PER SYMBOL
k = 1;

% SYMBOL RATE
rs = 28e9;

% NOISE BANDWIDTH 
bn = 12.5e9;

refbit = randi([0 1],k,nSymbol);

% mapping bit to symbol
sym = symbolizerGrayQam(refbit);

% symbol power
ps = sum(abs(sym).^2)/nSymbol; 

osnr = 0:16; % in dB
cr = 1; % coding rate
sps = 1;

for ndxOSNR = 1:length(osnr)
    
    snr(ndxOSNR) = osnr2snr(osnr(ndxOSNR),rs,sps,'complex');
    
    pn = ps/idbw(snr(ndxOSNR));
    
    signal = sym + genWGN(size(sym,1),size(sym,2),pn,'linear','complex');;

    bit = slicerGrayQam(signal,2^k);

    ber(ndxOSNR) = nnz(bit(:)-refbit(:))/(k*nSymbol);
    disp(sprintf('osnr = %d, ber = %.2e',osnr(ndxOSNR),ber(ndxOSNR)));
end
t_ber = T_BER_mQAM(osnr,2^k,rs);
% t_ber = (1/k)*3/2*erfc(sqrt(k*0.1*(10.^(EbNo_db/10))));
figure; plot(snr,log10(ber),'s-',snr,log10(t_ber),'k-'); grid on

