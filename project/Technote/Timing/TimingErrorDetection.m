%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This type of timing error detection uses a very simple Fourier transform
% property, i.e., time shift corresponds a phase ramp in spectrum.
% Therefore, the imaginary part of auto-correlation of spectrum with a
% fixed lag would be a TED with sinusoid characteristics. The smaller the
% lag, the longer range of detectable timing error.
%
% If the lag of auto-correlation is equal to the symbol rate, the range of
% detectable timing error is within one symbol duration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all


%% global parameters
% number of DFT
nfft = 4096;
% freq of signal
f1 = 1000;
f2 = 2000;
% sampling speed
fs = 8000;
%
nsample = 3850;


%% signal with 2 distinct frequencies
ntau = -4 : 0.1 : 4;
for ii = 1 : length(ntau)
    signal_1 = exp(1i * 2 * pi * f1 * ((0 : nsample - 1) + ntau(ii)) ./ fs);
    signal_2 = exp(1i * 2 * pi * f2 * ((0 : nsample - 1) + ntau(ii)) ./ fs);
    signal = signal_1 + signal_2;
    
    % get the periodogram, zero-padding the signal to nfft, power of 2
    signal_zp = [signal, zeros(1, nfft - nsample)];
    
    SS = fft(signal_zp) / nfft;
    
    ndx1 = 513;
    ndx2 = 1025;
    
    range = -15 : 15;
    timing_error(ii) = imag(sum(SS(range + ndx1) .* conj(SS(range + ndx2))));
end
figure; stem(ntau / (fs / f1), timing_error); grid on
xlabel('Timing error');
ylabel('TED output');

