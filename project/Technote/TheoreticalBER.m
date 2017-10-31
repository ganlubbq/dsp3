%% This script calculates the theoretical ber of m-qam signal in gaussian
% noise using one sample per symbol only, i.e., no pulse-shaping is
% considered in this version.

function [] = TheoreticalBER(input_k)
if nargin < 1, input_k = 1; end

ndata = 10^6;
% NUMBER OF BIT PER SYMBOL
k = input_k;
refbit = randi([0 1], k, ndata);
% mapping bit to symbol
sym = symbolizer_mqam(refbit);
% symbol power, assuming unit symbol period
ps = sum(abs(sym).^2) / ndata;

snr = -10 : 10; % in dB
for ndx = 1 : length(snr)
    pn = ps / idbw(snr(ndx));
    signal = sym + gaussian_noise(size(sym,1), size(sym,2), pn, 'linear', 'complex');
    bit = slicer_mqam(signal, 2^k);
    ber(ndx) = sum(abs(bit(:) - refbit(:))) / (k * ndata);
    fprintf('snr = %d, ber = %.2e\n', snr(ndx), ber(ndx));
end

if k == 1 
    % count only the noise in one dimension
    t_ber = snr2ber(idbw(snr) * 2, k, 'linear');
else
    t_ber = snr2ber(idbw(snr), k, 'linear');
end

figure; 
semilogy(snr, ber, 's-', snr, t_ber, 'k-'); grid on;
xlabel('SNR [dB]'); 
ylabel('BER'); 
legend(sprintf('%d bit per symbol', k), 'Theory');

return
