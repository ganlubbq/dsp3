% TEST SCRIPT FOR EVALUATING THE IMPLEMENTATION OF LINEAR CONVOLUTION BY
% USING FFT-BASED CIRCULAR CONVOLUTION. 
%
% Assuming that the length of filter l is always less than the length of
% data p, then padding both the filter coeff and data to l + p - 1 and do
% fft-based circular convolution would generate identical output with
% linear convolution for full length l + p - 1
%
% Oppenheim, Alan V., and Ronald W. Schafer. Discrete-time signal
% processing. Pearson Higher Education, 2010.
%
clear

% ***************************** Basic setting
% filter has length P < L
P = 6;
h = (1 : P) ./ P;

% data has length L
L = 11;
x = randi([0 1], 1, L);

% canonical conv output has length L+P-1
y = conv(x, h);

% circular conv will be identical to linear conv provided N >= L+P-1
N = L + P - 1;

% implement linear conv using N-point DFT (circular) and zero padding
yc_N = ifft(fft([x, zeros(1, N-L)]) .* fft([h, zeros(1, N-P)]));

% part of circular conv is identical to linear conv, if P < L, when using
% L-point DFT -- the first P-1 points of circular conv will be different
yc_L = ifft(fft(x) .* fft([h, zeros(1, L-P)]));

% to confirm...
h1 = figure; 
hold on; 
plot(y,'x-'); 
plot(1:N,yc_N,'gs-'); 
plot(yc_L,'rd-','MarkerFaceColor','r'); 
grid on; 
box on
hold off;
legend('Linear conv.','N-point DFT','L-point DFT');
title(sprintf('N = %d, L = %d, P = %d', N, L, P));

% ***************************** Long sequence - overlap and add
% filter impulse response, P=7
h = [1 2 3 4 3 2 1];
% data sequence, length 60
x = randi([0 1], 1, 60);
% linear convolution, length 66
y = conv(x, h);
% output
yy = zeros(1, 80);
% for each segment of data
L = 10;
% linear convolution result for each segment
N = 16;

% move L in each step and perform N-point circular convolution on L points
% of data based on zero padding and DFT, overlapping happens in output
for ii = 1 : 6
	xx = x((1 : L) + (ii - 1) * L);
	yy((1 : N) + (ii - 1) * L) = yy((1 : N) + (ii - 1) * L) + ifft(fft([xx, zeros(1,6)]) .* fft([h, zeros(1,9)]));
end

h2 = figure; 
hold on
plot(y,'x-'); 
plot(yy(1:length(y)),'o-'); 

% ***************************** Long sequence - overlap and save
% use the same data
yy = zeros(1,80);
L = 10;
N = 16;
p = 7;
x = [zeros(1, p-1), x];

% move L in each step and perform N-point circular convolution on N points
% of data based on DFT, and saving the last L points of output. Overlapping
% happens in input
for ii = 1 : 6
	xx = x((1 : N) + (ii - 1) * L);
	tmp = ifft(fft(xx) .* fft([h, zeros(1,9)]));
	yy((1 : L) + (ii - 1) * L) = tmp(p : end);
end

plot(yy(1:length(y)),'s-'); 
grid on
box on
hold off
legend('Linear conv.', 'overlap-add', 'overlap-save');

mngFigureWindow(h1, h2);
