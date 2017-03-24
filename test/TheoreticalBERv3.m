%% TEST SCRIPT FOR CALCULATING THEORETICAL BER VS OSNR
% CALCULATING THEORETICAL BER OF M-QAM SYSTEM WITH GAUSSIAN NOISE
% Denoted to pulse shaping, multiple samples per symbol, observing the ISI
% effet.

%%
function [] = TheoreticalBERv3(input_k)

if nargin < 1;
    input_k = 1;
end

nSymbol = 2^16;

% NUMBER OF BIT PER SYMBOL
k = input_k;

refbit = randi([0 1], k, nSymbol);

% mapping bit to symbol
sym = symbolizerGrayQam(refbit);

% symbol power
signal_power = sum(abs(sym).^2) / nSymbol; 

snr = -10 : 0.5 : 10; % in dB

sps = 16;

nSamples = sps * nSymbol;

% Upsampling
if ~iscolumn(sym)
    sym = sym.';
end
sym_upsampled = upSampInsertZeros(sym, sps);

% get a freq domain raised cosine filter response
Rs = 1;
Fs = sps;
% freqVect = getFFTGrid(nSamples,Fs);
alpha = 0.35; mode = 0;
H = calcRCFreqResponse(nSamples, Fs, Rs, alpha, mode);

% filtering signal in frequency domain
sym_upsampled_i = real(ifft(fft(real(sym_upsampled)) .* H));
sym_upsampled_q = real(ifft(fft(imag(sym_upsampled)) .* H));

sym_filtered = sym_upsampled_i + 1i * sym_upsampled_q;

signal_power_rcos_freq = sum(abs(sym_filtered).^2) / nSamples; % symbol power

% h1 = plotEyeDiagram(sym_filtered(1:8192),2*sps,'e');

data = sym_filtered(1:sps:end);
symbol_power = sum(abs(data).^2) / nSymbol;

snr = -10 : 10; % in dB

for ndx = 1:length(snr)
    
    % noise power
    pn = symbol_power / idbw(snr(ndx));
    
    % awgn
    signal = data + genWGN(size(data,1), size(data,2), pn, 'linear', 'complex');
   
    bit = slicerGrayQam(signal, 2^k);
   
    ber(ndx) = sum(abs(bit(:) - refbit(:))) / (k*nSymbol);
    
    disp(sprintf('snr = %d, ber = %.2e', snr(ndx), ber(ndx)));
end

if k == 1 
    % count only real noise 
    t_ber = T_BER_SNR_mQAM(idbw(snr)*2, 2^k);
else
    t_ber = T_BER_SNR_mQAM(idbw(snr), 2^k);
end

figure; grid on
plot(snr,log10(ber),'s-',snr,log10(t_ber),'k-'); 
xlabel('SNR dB'); ylabel('LOG10 BER'); 
legend(sprintf('%d bit per symbol',k),'Theory');

return

