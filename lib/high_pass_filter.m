% DESCRIPTION
% 
% Example: 
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright 2015 default
function y = high_pass_filter(x,Rs,Fs)

beta = 0.9999

T = 1 / Rs;

N = length(x);
n = 1 : N/2;
RF = Fs / (N-1) / 2 * (2*n-1);
LF = fliplr(RF) * -1.0;
f = [LF,RF];

H = zeros(size(f));

for ii = 1:length(f)
    if abs(f(ii)) <= (1-beta)/(2*T)
        H(ii) = T;
    elseif (abs(f(ii))>(1-beta)/(2*T) && abs(f(ii))<=(1+beta)/(2*T))
        H(ii) = (T/2) * (1+cos(pi*T/beta*(abs(f(ii))-(1-beta)/(2*T))));
    else
        H(ii) = 0;
    end
end

H = H(:)/T;

detune_ratio = 0.15;

delta_f = detune_ratio * Rs;

delta_f_N = delta_f / Fs * N;

H = circshift(H,delta_f_N);

F = fftshift(fft(x));

FH = F.*H;

y = ifft(ifftshift(FH));

% figure; plot(log10(abs(FH).^2))

return
