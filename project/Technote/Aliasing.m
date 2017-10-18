clear
close all

nsample = 3999;
fs = 200;
order = 0.35;
bandwidth = 2;
H = frequency_response(nsample, fs, order, bandwidth, 'rc');
freq = getFFTGrid(nsample, fs);
% figure; plot(fftshift(freq), fftshift(H), 'k-'); grid on; hold on; box on;

% simulate a continous-time signal x(t) bandlimited to -0.5 ~ 0.5
h = ifftshift(ifft(H));
figure; plot(1:nsample, real(h)); grid on; box on; hold on


% sampling speed 1.35 * 2 - nyquist
ac_fs = 1.35 * 2;
hs = h(1 : fs/ac_fs : end);
stem(1 : fs/ac_fs : nsample, hs);

% low pass filter with bw of 1.35
H1 = frequency_response(nsample, fs, [], 1.35, 'rect');
hs = upSampInsertZeros(hs, fs/ac_fs);
hr = ifft(fft(hs(1:nsample)) .* H1);
hr = hr .* (max(h) / max(hr));
plot(1:nsample, real(hr)); grid on; box on 



% sampling speed 4 - oversampling
figure; plot(1:nsample, real(h)); grid on; box on; hold on
ac_fs = 4;
hs = h(1 : fs/ac_fs : end);
stem(1 : fs/ac_fs : nsample, hs);

% low pass filter with bw of 0.5 ~ 2
H1 = frequency_response(nsample, fs, [], 0.7, 'rect');
hs = upSampInsertZeros(hs, fs/ac_fs);
hr = ifft(fft(hs(1:nsample)) .* H1);
hr = hr .* (max(h) / max(hr));
plot(1:nsample, real(hr)); grid on; box on 



% sampling speed 0.5 - undersampling
figure; plot(1:nsample, real(h)); grid on; box on; hold on
ac_fs = 0.4;
hs = h(1 : fs/ac_fs : end);
stem(1 : fs/ac_fs : nsample, hs);

% low pass filter with bw of 0.25
H1 = frequency_response(nsample, fs, [], 0.25, 'rect');
hs = upSampInsertZeros(hs, fs/ac_fs);
hr = ifft(fft(hs(1:nsample)) .* H1);
hr = hr .* (max(h) / max(hr));
plot(1:nsample, real(hr)); grid on; box on 










