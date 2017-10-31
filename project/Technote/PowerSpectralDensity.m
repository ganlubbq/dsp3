% Power spectral density of linear digital modulation contains two parts:
% (1) the power spectral of pulse-shaping filter; (2) discrete harmonics
% proportional to the mean of information.
%
% Ref: Digital Communications By John Proakis 5th Edition
clear;
% ASK
nbaud = 8000; % sample size
H = frequency_response(nbaud * 16, 16, 3, .9, 'bessel');
for ii = 1:20
    baud = randi([0, 1], nbaud, 1);
    baud_wfm = real(ifft(fft(upsampling(baud, 16, 1)) .* H));
    psd(:,ii) = abs(fft(baud_wfm)).^2 / numel(baud_wfm);
end
psd = mean(psd, 2);
figure; plot(baud_wfm(1:160)); grid on;
figure; 
semilogy(getFFTGridPos(nbaud * 16, 16), psd(1:end/2)); hold on;
semilogy(getFFTGridPos(nbaud * 16, 16), abs(H(1:end/2)).^2/100, 'LineWidth', 4); grid on;
xlabel('Frequency [Hz]');
ylabel('Power spectral density [dB]');
% PSK
nbaud = 8000; % sample size
H = frequency_response(nbaud * 16, 16, 3, .9, 'bessel');
for ii = 1:20
    baud = randi([0, 1], nbaud, 1) - 0.5;
    baud_wfm = real(ifft(fft(upsampling(baud, 16, 1)) .* H));
    psd(:,ii) = abs(fft(baud_wfm)).^2 / numel(baud_wfm);
end
psd = mean(psd, 2);
figure; plot(baud_wfm(1:160)); grid on;
figure; 
semilogy(getFFTGridPos(nbaud * 16, 16), psd(1:end/2)); hold on;
semilogy(getFFTGridPos(nbaud * 16, 16), abs(H(1:end/2)).^2/100, 'LineWidth', 4); grid on;
xlabel('Frequency [Hz]');
ylabel('Power spectral density [dB]');
% QAM