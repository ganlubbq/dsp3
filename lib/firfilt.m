% Implementing FIR filter using overlap-and-add and overlap-and-save
% structures
% 
% Usage: y = firfilt(x,h,method,shape) filters input data x with filter
% impulse response in time domain h. The implementation method could either
% be 'overlap-add' or 'overlap-save'
% 
% Input: 
%       x       - input data with length L
%       h       - filter coeff. with length P < L
%       y       - output data with max length of L+P-1
%       method  - structure of FIR filter
%       shape   - shape of output data
% 
% Reference: 
% 
% Note: The input and output are in row format. Length of filter coeff is
%   less than length of input data. FFT length is at least 2x filter length
% 
% See Also: conv, filter, fftfilt.
%
% Copyright 2015 DAWEI DE

function y = firfilt(x,h,method,shape)

if nargin<4
    shape = 'full';
end

if nargin<3
    method = 'overlap-add';
end

if ~isrow(x) || ~isrow(h)
    error('The first and second inputs have to be row vectors');
end


% default convolution output with length of L+P-1
yc = conv(x,h);

switch method
    
    case 'overlap-add'
        
        % divide x into sectoins with length of L
        % h has length of P
        % FFT is carried out with length of L+P-1 which is power of 2
        
        P = length(h);
        
        N = 2^nextpow2(P*2);
        
        L = N+1-P;
        
        cutoff = L - mod(length(x),L);
        
        % pending 0s to x to be multiple length of L
        x = [x,zeros(1,cutoff)];
        
        yc_N = zeros(1,length(x)+P-1);
        
        for ii = 1:length(x)/L
            % move with step of L
            xx = x((1:L)+(ii-1)*L);
            % implementing overlap and add method
            yc_N((1:N)+(ii-1)*L) = yc_N((1:N)+(ii-1)*L) + ifft(fft([xx zeros(1,P-1)]).*fft([h zeros(1,L-1)]));
        end
        
        yy = yc_N(1:end-cutoff);
        
    case 'overlap-save'
        P = length(h);
        
        N = 2^nextpow2(P*2);
        
        L = N;
        
        yc_L = zeros(1,length(x)+P-1);
        
        s = L-P+1;

        x = [zeros(1,P-1) x];
        
        cutoff = s - mod(length(x),s);
        % pending 0s to x to be multiple length of s: n*s+L
        x = [x,zeros(1,cutoff),zeros(1,P-1)];
        
        for ii = 1:floor(length(x)/s)
            xx = x((1:L)+(ii-1)*s);
            tmp = ifft(fft(xx).*fft([h zeros(1,L-P)]));
            % tmp = ifft(fft(xx).*fft(h));
            yc_L((1:s)+(ii-1)*s) = tmp(P:end);
        end
        yy = yc_L(1:end-cutoff);
        
    otherwise
        error('unsupported FIR filter method!!!');
end

switch shape
    case 'full'
        y = yy;
    case 'same'
        y = yy(1:end-P+1);
    case 'valid'
        y = yy(P:end-P+1);
    otherwise
        error('unsupported FIR filter shape!!!');
end

% figure; hold on; grid on; 
% plot(yc,'o-'); 
% plot(y,'gs-'); 
% legend('linear convolution',sprintf('%s',method)); hold off;

return