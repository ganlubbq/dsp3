%% TEST SCRIPT FOR CALCULATING THEORETICAL BER VS OSNR
% CALCULATING THEORETICAL BER OF M-QAM SYSTEM WITH GAUSSIAN NOISE
%%
clear

nSymbol = 10^6;

% NUMBER OF BIT PER SYMBOL
k = 1

refbit = randi([0 1],k,nSymbol);

% mapping bit to symbol
sym = symbolizerGrayQam(refbit);

% symbol power
ps = sum(abs(sym).^2)/nSymbol; 

snr = -10:10; % in dB
cr = 1; % coding rate
sps = 1;

for ndx = 1:length(snr)
    
    % noise power
    pn = ps/idbw(snr(ndx));
    
    % awgn
    signal = sym + genWGN(size(sym,1),size(sym,2),pn,'linear','complex');
   
    bit = slicerGrayQam(signal,2^k);
   
    ber(ndx) = nnz(bit(:)-refbit(:))/(k*nSymbol);
    
    disp(sprintf('snr = %d, ber = %.2e',snr(ndx),ber(ndx)));
end

if k==1 % count only real noise 
    t_ber = T_BER_SNR_mQAM(idbw(snr)*2,2^k);
else
    t_ber = T_BER_SNR_mQAM(idbw(snr),2^k);
end

figure; plot(snr,log10(ber),'s-',snr,log10(t_ber),'k-'); grid on
xlabel('SNR dB'); ylabel('LOG10 BER');

