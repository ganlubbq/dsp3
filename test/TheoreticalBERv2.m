% Test script for calculating theoretical ber of m-qam system with gaussian
% noise, using one sample per symbol only, i.e. no pulse shaping in this
% version

% Using input argument to specify the number of bits per symbol
function [] = TheoreticalBERv2(input_k)
if nargin < 1;
    input_k = 1;
end

nSymbol = 10^6;

% NUMBER OF BIT PER SYMBOL
k = input_k;

refbit = randi([0 1], k, nSymbol);

% mapping bit to symbol
sym = symbolizerGrayQam(refbit);

% symbol power, assuming unit symbol period
ps = sum(abs(sym).^2) / nSymbol; 

snr = -10 : 10; % in dB
for ndx = 1 : length(snr)
    % noise power
    pn = ps / idbw(snr(ndx));
    
    % awgn
    signal = sym + genWGN(size(sym,1), size(sym,2), pn, 'linear', 'complex');
   
    bit = slicerGrayQam(signal, 2^k);
   
    ber(ndx) = sum(abs(bit(:) - refbit(:))) / (k * nSymbol);
    
    fprintf('snr = %d, ber = %.2e\n', snr(ndx), ber(ndx));
end

if k == 1 
    % count only the noise in one dimension
    t_ber = snr2ber(idbw(snr) * 2, k, 'linear');
else
    t_ber = snr2ber(idbw(snr), k, 'linear');
end

figure; 
semilogy(snr, ber, 's-', snr, t_ber, 'k-');
grid on;
xlabel('SNR [dB]'); 
ylabel('BER'); 
legend(sprintf('%d bit per symbol', k), 'Theory');

return
