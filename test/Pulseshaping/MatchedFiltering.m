% Matched filtering is the max-snr filtering technique. The signal model
%       r = x + n
% the filter shape acting on r that could get the maximum snr is exactly
% the x so that the signal is recovered as the self inner product.

clear

nSymbol = 2^10;

% NUMBER OF BIT PER SYMBOL
k = 2;
refbit = randi([0 1], k, nSymbol);

% mapping bit to symbol
sym = symbolizerGrayQam(refbit);

% symbol power
signal_power = sum(abs(sym).^2) / nSymbol

sps = 16;
nSamples = sps * nSymbol;

%% Upsampling
if ~iscolumn(sym)
    sym = sym.';
end
sym_upsampled = sps * upSampInsertZeros(sym, sps);

%% get a freq domain raised cosine filter response
Rs = 1;
Fs = sps;
freqVect = getFFTGrid(nSamples,Fs);
alpha = 0.35; 
% RRC
mode = 1;
H = calcRCFreqResponse(nSamples, Fs, Rs, alpha, mode);

% filtering signal in frequency domain
sym_upsampled_i = real(ifft(fft(real(sym_upsampled)) .* H));
sym_upsampled_q = real(ifft(fft(imag(sym_upsampled)) .* H));

txSignal = sym_upsampled_i + 1i*sym_upsampled_q;

% display signal power after pulse shaping, optimal sampling instance
signal_power_tx = sum(abs(txSignal(1:sps:end)).^2) / nSymbol

h1 = plotEyeDiagram(txSignal(1:end), 2*sps, 'e');

% set a fixed snr dB
snr = 10;

noise_power = 2 / idbw(snr);
noise = genWGN(size(txSignal,1), size(txSignal,2), noise_power, 'linear', 'complex');

rxSignal = txSignal + noise;

% filtering received signal in frequency domain by using a matched filter
rxSignal_i = real(ifft(fft(real(rxSignal)) .* H));
rxSignal_q = real(ifft(fft(imag(rxSignal)) .* H));

recSignal = rxSignal_i + 1i * rxSignal_q;

h2 = plotEyeDiagram(recSignal(1:end), 2*sps, 'e');

mngFigureWindow(h1,h2);

% read the snr of recovered signal
recPS = sum(abs(recSignal(1:sps:end)).^2) / nSymbol;
% recSNR = dbw();
