function y = ted_FLN(px)
%FLN Fourth-Law Nonlinearity

N = length(px);

k = 1:N;
ex = exp(-1j.*(k-1).*pi./2);

s = sum( abs(px).^4 .* ex.' );

y = -angle(s) / 2 / pi;