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
close all

%% global parameter...
sps     = 60;
nsymbol = 1024;
nsample = sps * nsymbol;

k = 2;
refbit = randi([0 1], k, nsymbol);
symseq = symbolizerGrayQam(refbit);
% symseq = ones(1,nSymbol) - 1i * ones(1,nSymbol);

% energy and power of symbol sequence
T = 1;
symSeqEnergy = T * sum(abs(symseq).^2);
symSeqPower  = T * sum(abs(symseq).^2) / (T * nsymbol);

% energy and power of zero-padded symbol sequence (the signal)
Ts = T / sps;
pwr_sig = Ts * sum(abs(symseq).^2) / (Ts * nsample);

% By the definition, the energy of discrete signal is approximation of
% energy of its continuous version, i.e, replacing the integral of power
% with summation.

% The energy of symbol sequence goes smaller and smaller as the sampling
% speed increases. It can be understood as more and more "emptiness" is
% discovered with increasing sampling speed.

% Imaging that somehow we are able to sample the symbol sequence at the
% same speed as its repetition without aliasing, and incredibly we're
% sampling at the correct instance, we have no choice but to believe that
% the signal being sampled evolves slowly and smoothly. However, as the
% sampling speed goes faster, we start to realize that the signal being
% sampled is nothing but only a constant power occuring in a very short
% period of time, i.e., the sampling interval. Therefore, a smaller signal
% energy is concluded.

% Luckily, the zero-padded symbol sequence is an extremely non-practical in
% the real world as it has infinite bandwidth and there is no way to sample
% it properly. Signals in the real world are inevitably band-limited and
% therefore exhibit smoothness in time domain to some extent. Hence, a
% smaller sampling interval will give a better approximation of energy.


%% pulse-shaping in frequency domain...
beta = 0.65;

% upsampling
sym_upsampled = upSampInsertZeros(symseq(:), sps);

% freq response of raised cosine filter with unit energy, i.e.,
% sum(ifft(H).^2) == 1
% such that the signal after pulse-shaping will have the same average power
H = calcFilterFreqResp(nsample, sps, beta, 1, 'rc');
H = H./sqrt(sum(H.^2)/nsample);

% filtering signal in frequency domain
sym_upsampled_i = real(ifft(fft(real(sym_upsampled)) .* H));
sym_upsampled_q = real(ifft(fft(imag(sym_upsampled)) .* H));

sym_filtered = sym_upsampled_i + 1i * sym_upsampled_q;

% signal power after pulse-shaping, which is NOT equal to the power of
% optimal sampling instances, i.e., pwr_sig_freq' =
% calcrms(sym_filtered(1:sps:end)).^2. This power will be larger than the
% averaged power especially for large roll-off
pwr_sig_freq = calcrms(sym_filtered).^2;

h1 = plotEyeDiagram(sym_filtered(1:end), 2*sps, 'e');
title('Frequency domain filtering');


%% pulse-shaping in time domain...
beta = 0.65;

% impulse response of rcos digital filter with unit energy, i.e.,
% sum(h.^2) == 1
% such that the signal after pulse-shaping will have the same average power
span = 10;
h = rcosdesign(beta, span, sps, 'normal');

sym_upsampled_i = upfirdn(real(symseq(:)), h, sps);
sym_upsampled_q = upfirdn(imag(symseq(:)), h, sps);

% sym_upsampled_i = firfilt(real(sym_upsampled(:)), h, 'overlap-save','same');
% sym_upsampled_q = firfilt(imag(sym_upsampled(:)), h, 'overlap-save','same');

delay = span * sps / 2;

sym_filtered = sym_upsampled_i + 1i * sym_upsampled_q;

% signal power after pulse-shaping
pwr_sig_time = calcrms(sym_filtered).^2;

h2 = plotEyeDiagram(sym_filtered(1 + delay : end - delay), 2*sps, 'e');
title('Time domain filtering');

mngFigureWindow(h1, h2);


%% print results...
% The energy of one symbol is spread into one pulse by pulseshaping, if the
% pulseshaping filter has unit energy
fprintf('True signal power is %.2f dB\n', dbw(pwr_sig));
fprintf('Freq-rcos signal power is %.2f dB\n', dbw(pwr_sig_freq));
fprintf('Time-rcos signal power is %.2f dB\n', dbw(pwr_sig_time));

