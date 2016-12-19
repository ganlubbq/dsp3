function y = ted_SLN(px)
%SLN Square-Law Nonlinearity

N = length(px);

k = 1:N;
ex = exp(-1j.*(k-1).*pi./2);

s = sum( abs(px).^2 .* ex.' );

y = -angle(s) / 2 / pi;