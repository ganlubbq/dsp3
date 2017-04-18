% TEST SCRIPT FOR PULSE-SHAPING USING RCOS FILTER Compare raised cosine
% pulse-shaping filter in both frequency and time domain. Note that in both
% cases, the upsampling should be a zero inserting one.
clear

nSymbol = 2^10;

% NUMBER OF BIT PER SYMBOL
k = 2;

refbit = randi([0 1], k, nSymbol);

% mapping bit to symbol
sym = symbolizerGrayQam(refbit);
% sym = ones(1,nSymbol) - 1i * ones(1,nSymbol);

% symbol power
pwr_sig = sum(abs(sym).^2) / nSymbol;

sps = 64;

nSamples = sps * nSymbol;

% filter response domain, should be FREQ or TIME
domain = 'FREQ';
% domain = 'TIME';

% roll-off factor of raised cosine filter
alpha = 0.35;

%%%%%%%%%%%%%%%%%%%%%%% FREQ %%%%%%%%%%%%%%%%%%%%%%%
% upsampling
sym_upsampled = upSampInsertZeros(sym(:), sps);

% freq response of raised cosine filter with unit energy
H = calcRCFreqResponse(nSamples, sps, 1, alpha, 'rc');
H = H./sqrt(sum(H.^2)/nSamples);

% filtering signal in frequency domain
sym_upsampled_i = real(ifft(fft(real(sym_upsampled)) .* H));
sym_upsampled_q = real(ifft(fft(imag(sym_upsampled)) .* H));

delay = 0;

sym_filtered = sym_upsampled_i + 1i * sym_upsampled_q;

% display signal power after pulse shaping
pwr_sig_freq = sum(abs(sym_filtered).^2) / nSamples;
h1 = plotEyeDiagram(sym_filtered(1:end), 2*sps, 'e');

%%%%%%%%%%%%%%%%%%%%%%% TIME %%%%%%%%%%%%%%%%%%%%%%%
% impulse response of rcos digital filter
span = 10;
h = rcosdesign(alpha, span, sps, 'normal');

sym_upsampled_i = upfirdn(real(sym(:)), h, sps);
sym_upsampled_q = upfirdn(imag(sym(:)), h, sps);

% sym_upsampled_i = firfilt(real(sym_upsampled(:)), h, 'overlap-save','same');
% sym_upsampled_q = firfilt(imag(sym_upsampled(:)), h, 'overlap-save','same');

delay = span * sps / 2;

sym_filtered = sym_upsampled_i + 1i * sym_upsampled_q;

% display signal power after pulse shaping
pwr_sig_time = sum(abs(sym_filtered).^2) / nSamples;


sym_filtered = sym_filtered(:);
h2 = plotEyeDiagram(sym_filtered(1 + delay : end - delay), 2*sps, 'e');

mngFigureWindow(h1, h2);

% The energy of one symbol is spread into one pulse by pulseshaping, if the
% pulseshaping filter has unit energy
fprintf('True signal power is %.2f dB\n', dbw(pwr_sig/sps));
fprintf('Freq-rcos signal power is %.2f dB\n', dbw(pwr_sig_freq));
fprintf('Time-rcos signal power is %.2f dB\n', dbw(pwr_sig_time));

