% This script calculates the theoretical ber of m-qam signal in gaussian
% noise using multiple sample per symbol with pulse-shaping considered in
% this version. It is shown that the real SNR is obtained by doing matched
% filtering at the receiver.
function [] = MatchedFiltering(input_k)
if nargin < 1;
    input_k = 1;
end

sps = 16;

nSymbols = 2^16;
nSamples = sps * nSymbols;

% NUMBER OF BIT PER SYMBOL
k = input_k;

refbit = randi([0 1], k, nSymbols);

% mapping bit to symbol
sym = symbolizerGrayQam(refbit);

% pulse-shaping using frequency domain method
% get a freq domain raised cosine filter response
alpha = 0.35;
H = calcFilterFreqResp(nSamples, sps, alpha, 1, 'rrc');

sym_upsampled = upSampInsertZeros(sym(:), sps);

% filtering signal in frequency domain
sym_upsampled_i = real(ifft(fft(real(sym_upsampled)) .* H));
sym_upsampled_q = real(ifft(fft(imag(sym_upsampled)) .* H));

txSignal = sym_upsampled_i + 1i * sym_upsampled_q;

% get the signal power which is the power of information
symbol_power = calcrms(txSignal).^2;

% this is symbol snr in db
snr = -10 : 1 : 10;
for ndx = 1 : length(snr)
    % noise power in the whole spectrum, obtained by multiplying the noise
    % power per symbol by the oversampling factor
    pn = symbol_power * sps / idbw(snr(ndx));
    
    % awgn
    rxSignal = txSignal + genWGN(size(txSignal,1), size(txSignal,2), pn, 'linear', 'complex');

    % matched filtering
    signal_i = real(ifft(fft(real(rxSignal)) .* H));
    signal_q = real(ifft(fft(imag(rxSignal)) .* H));
    rxSym = signal_i(1:sps:end) + 1i * signal_q(1:sps:end);
    
    bit = slicerGrayQam(rxSym, 2^k);
   
    ber(ndx) = sum(abs(bit(:) - refbit(:))) / numel(refbit);
    
    fprintf('snr = %.1f, ber = %.2e\n', snr(ndx), ber(ndx));
end

if k == 1 
    % count only real noise 
    t_ber = snr2ber(idbw(snr)*2, k, 'linear');
else
    t_ber = snr2ber(idbw(snr), k, 'linear');
end

figure;
semilogy(snr, ber, 's-', snr, t_ber, 'k-'); grid on
xlabel('SNR dB'); ylabel('LOG10 BER'); 
legend(sprintf('%d bit per symbol',k),'Theory');

return
