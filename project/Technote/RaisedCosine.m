clear
close all

x = -6 : 0.01 : 6;
y = sinc(x);

figure; plot(x,y); grid on; box on; hold on
x = -6 : 1 :6;
y = sinc(x);
stem(x,y); grid on; box on


x = -6 : 0.01 : 6;
y = sinc(x);
figure; plot(x, y, 'k', 'LineWidth', 1); grid on; box on; hold on

alpha = 0.202;
T = 1;
t = -6 : 0.01 : 6;
h = sin(pi*t./T) ./ (pi*t/T) .* cos(alpha * pi * t / T) ./ (1 - 4 * alpha^2 * t.^2 / T^2);
plot(t, h, 'k--', 'LineWidth', 1);

alpha = 0.35;
T = 1;
t = -6 : 0.01 : 6;
h = sin(pi*t./T) ./ (pi*t/T) .* cos(alpha * pi * t / T) ./ (1 - 4 * alpha^2 * t.^2 / T^2);
plot(t, h, 'k-.', 'LineWidth', 1);

x = -6 : 6;
y = sinc(x);
stem(x, y, 'k', 'LineWidth', 2);

legend('Perfect SINC', 'Raised Cosine 0.2', 'Raised Cosine 0.35', 'Samples');


nsample = 2^12;
fs = 2;
order = 0.00002;
bandwidth = 1;
type = 'rc';
H = frequency_response(nsample, fs, order, bandwidth, type);
freq = getFFTGrid(nsample, fs);
figure; plot(fftshift(freq), fftshift(H), 'k-'); grid on; hold on; box on;

order = 0.2;
bandwidth = 1;
type = 'rc';
H = frequency_response(nsample, fs, order, bandwidth, type);
freq = getFFTGrid(nsample, fs);
plot(fftshift(freq), fftshift(H), 'k--'); 

order = 0.35;
bandwidth = 1;
type = 'rc';
H = frequency_response(nsample, fs, order, bandwidth, type);
freq = getFFTGrid(nsample, fs);
plot(fftshift(freq), fftshift(H), 'k-.'); 

ylim([-0.2, 1.2]);

legend('Perfect Rectangular', 'Raised Cosine 0.2', 'Raised Cosine 0.35');
