function y = ted_LOGN(px)
%LOGN LOG Nonlinearity

N = length(px);

k = 1:N;
ex = exp(-1j.*(k-1).*pi./2);

f = log(1+10.*abs(px).^2);

s = sum( f .* ex.' );

y = -angle(s) / 2 / pi;