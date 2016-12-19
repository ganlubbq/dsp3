function y = ted_lee(x,g)
% Yan Wang, "An Alternative Blind Feddforward Symbol Timing Estimator Using
% Two Samples per Symbol," IEEE TRANSACTIONS ON COMMUNICATIONS, VOL. 51,
% NO. 9, SEPTEMBER 2003

if nargin<2
    g = 1.414;
end

L = length(x);

n = 1:L;

% ex1 = exp(-1j.*(ii-1).*pi);
ex1 = (-1).^(n-1);

sum_1 = sum( abs(x).^2 .* ex1.' );

% ex2 = exp(-1j.*(ii-1.5).*pi);
ex2 = 1j * (-1).^(n-1);

xh = x(2:end);
xx = x(1:end-1);
ex2 = ex2(1:end-1);

sum_2 = sum( real(conj(xx).*xh) .* ex2.' );

sum_3 = g*sum_1 + sum_2;

y = angle(sum_3) / 2 / pi;

