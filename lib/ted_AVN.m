function y = ted_AVN(px)

N = length(px);

k = 1:N;
ex = exp(-1j.*(k-1).*pi./2);

s = sum( abs(px) .* ex.' );

y = -angle(s) / 2 / pi;