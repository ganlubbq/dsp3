function y = ted_SLN_int(px)
%SLN Square-Law Nonlinearity

N = length(px);

x = interpft(px, 2*N);

k = 1:2*N;
ex = exp(-1j.*(k-1).*pi./2);

s = sum( abs(x).^2 .* ex.' );

y = -angle(s) / 2 / pi;