function y = ted_gar(px)
%GARDNER Gardner timing error detector
%
% The modified version is refered to:
%
% [1] W. Gappmair, S. Cioni, G. E. Corazza, and O. Koudelka, ¡°Symbol-Timing
% Recovery with Modified Gardner Detectors,¡± in International Symposium on
% Wireless Communication Systems, 2005, pp. 831¨C834.
%

px1= px(1:2:end-2);
px2= px(2:2:end-1);
px3= px(3:2:end);

y = mean( (real(px1)-real(px3)).*real(px2) + ...
           (imag(px1)-imag(px3)).*imag(px2) );

% mu = -0.5;
% y  = mean( (real(px3).* (abs(px3).^ (mu-1))-real(px1).* (abs(px1).^ (mu-1))).*real(px2) + ...
%            (imag(px3).* (abs(px3).^ (mu-1))-imag(px1).* (abs(px1).^ (mu-1))).*imag(px2) );

% px3 = abs(px3).^ mu.* exp(1i*angle(px3));
% px1 = abs(px1).^ mu.* exp(1i*angle(px1));
% y = mean( real((px3-px1).* conj(px2)) );