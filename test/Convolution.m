%% TEST SCRIPT FOR EVALUATING THE IMPLEMENTATION OF LINEAR CONVOLUTION BY USING CIRCULAR CONVOLUTION BASED ON DFT
% ASSUMING LENGTH OF FILTER L IS ALWAYS LESS THAN THAT OF DATA P, PADDING
% FILTER COEFF AND DATA BOTH TO L+P-1 AND DO DFT-BASED CIRCULAR CONVOLUTION
% WOULD GENERATE IDENTICAL OUTPUT WITH LINEAR CONVOLUTION FOR FULL LENGTH
%
% Oppenheim, Alan V., and Ronald W. Schafer. Discrete-time signal
% processing. Pearson Higher Education, 2010.
%
%%
clear
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',148));

%% Basic setting

% filter has length P
P = 6;
h = (1:P)./P;

% data has length L
L = 11;
x = randi([0 1],1,L);

% canonical conv output has length L+P-1
y = conv(x,h);

% circular conv will be identical to linear conv provided N >= L+P-1
N = L + P - 1;

% implement linear conv using N-point DFT (circular) and zero padding
yc_N = ifft(fft([x zeros(1,N-L)]) .* fft([h zeros(1,N-P)]));

% part of circular conv is identical to linear conv if P < L when using
% L-point DFT -- the first P-1 points of circular conv will be different
yc_L = ifft(fft(x) .* fft([h zeros(1,L-P)]));

% to confirm...
figure; title(sprintf('N = %d, L = %d, P = %d',N,L,P));
hold on; 
plot(y,'x-'); 
plot(1:N,yc_N,'gs-'); 
plot(yc_L,'rd-','MarkerFaceColor','r'); grid on; legend('conv','N-point DFT','L-point DFT');
hold off;

%% Long sequence - overlap and add
% filter impulse response, P=7
h = [1 2 3 4 3 2 1];
% data sequence, length 60
x = randi([0 1],1,60);
% linear convolution, length 66
y = conv(x,h);
% output
yy = zeros(1,80);
% for each segment of data
L = 10;
% linear convolution result for each segment
N = 16;

% move L in each step and perform N-point circular convolution on L points
% of data based on zero padding and DFT, overlapping happens in output
for ii = 1:6
	xx = x((1:L) + (ii-1) * L);
	yy((1:N) + (ii-1) * L) = yy((1:N) + (ii-1) * L) + ifft(fft([xx zeros(1,6)]) .* fft([h zeros(1,9)]));
end

figure; hold on
plot(y,'x-'); plot(yy(1:length(y)),'o-'); grid on

%% Long sequence - overlap and save
% filter P=7
h = [1 2 3 4 3 2 1];
% data sequence
x = randi([0 1],1,60);
% linear convolution with length of 66
y = conv(x,h);

yy = zeros(1,80);

L = 10;

N = 16;

p = 7;

x = [zeros(1,p-1) x];

% move L in each step and perform N-point circular convolution on N points
% of data based on DFT, and saving the last L points of output. Overlapping
% happens in input
for ii = 1:6
	xx = x((1:N) + (ii-1) * L);
	tmp = ifft(fft(xx) .* fft([h zeros(1,9)]));
% tmp = ifft(fft(xx).*fft(h));
	yy((1:L) + (ii-1) * L) = tmp(p:end);
end

figure; hold on
plot(y,'x-'); plot(yy(1:length(y)),'o-'); grid on

