% The raised cosine pulse-shaping filters in the frequency and time domain
% are compared. In both cases, the type of upsampling symbol sequence
% should be zero inserting. It is shown that the pulse-shaping filter with
% unit energy (integral of instantaneous power) preserves the absolute
% power of signal, i.e., the power of information. The RC filter with small
% roll-off factor has large PAPR such that downsampling the signal to 1 sps
% with different offset will result in the same power. On the other hand,
% the RC filter with large roll-off has small PAPR such that the power of
% downsampled signal will be maximized at the optimal sampling instances.
% All results apply to RRC filter as well.
clear

sps = 64;

nSymbols = 2^10;
nSamples = sps * nSymbols;

% NUMBER OF BIT PER SYMBOL
k = 2;

refbit = randi([0 1], k, nSymbols);

% mapping bit to symbol
sym = symbolizerGrayQam(refbit);
% sym = ones(1,nSymbol) - 1i * ones(1,nSymbol);

% symbol power spread into one symbol period after zero inserting
pwr_sig = calcrms(sym).^2 / sps;

% filter response domain, should be FREQ or TIME
domain = 'FREQ';
% domain = 'TIME';

% roll-off factor of raised cosine filter
beta = 0.65;

%%%%%%%%%%%%%%%%%%%%%%% FREQ 
% upsampling
sym_upsampled = upSampInsertZeros(sym(:), sps);

% freq response of raised cosine filter with unit energy, i.e.,
% sum(ifft(H).^2) == 1
% such that the signal after pulse-shaping will have the same average power
H = calcRCFreqResponse(nSamples, sps, 1, beta, 'rc');
H = H./sqrt(sum(H.^2)/nSamples);

% filtering signal in frequency domain
sym_upsampled_i = real(ifft(fft(real(sym_upsampled)) .* H));
sym_upsampled_q = real(ifft(fft(imag(sym_upsampled)) .* H));

delay = 0;

sym_filtered = sym_upsampled_i + 1i * sym_upsampled_q;

% signal power after pulse-shaping, which is NOT equal to the power of
% optimal sampling instances, i.e., pwr_sig_freq' =
% calcrms(sym_filtered(1:sps:end)).^2. This power will be larger than the
% averaged power especially for large roll-off
pwr_sig_freq = calcrms(sym_filtered).^2;

h1 = plotEyeDiagram(sym_filtered(1:end), 2*sps, 'e');

%%%%%%%%%%%%%%%%%%%%%%% TIME 
% impulse response of rcos digital filter with unit energy, i.e.,
% sum(h.^2) == 1
% such that the signal after pulse-shaping will have the same average power
span = 10;
h = rcosdesign(beta, span, sps, 'normal');

sym_upsampled_i = upfirdn(real(sym(:)), h, sps);
sym_upsampled_q = upfirdn(imag(sym(:)), h, sps);

% sym_upsampled_i = firfilt(real(sym_upsampled(:)), h, 'overlap-save','same');
% sym_upsampled_q = firfilt(imag(sym_upsampled(:)), h, 'overlap-save','same');

delay = span * sps / 2;

sym_filtered = sym_upsampled_i + 1i * sym_upsampled_q;

% signal power after pulse-shaping
pwr_sig_time = calcrms(sym_filtered).^2;

h2 = plotEyeDiagram(sym_filtered(1 + delay : end - delay), 2*sps, 'e');

mngFigureWindow(h1, h2);

% The energy of one symbol is spread into one pulse by pulseshaping, if the
% pulseshaping filter has unit energy
fprintf('True signal power is %.2f dB\n', dbw(pwr_sig));
fprintf('Freq-rcos signal power is %.2f dB\n', dbw(pwr_sig_freq));
fprintf('Time-rcos signal power is %.2f dB\n', dbw(pwr_sig_time));

