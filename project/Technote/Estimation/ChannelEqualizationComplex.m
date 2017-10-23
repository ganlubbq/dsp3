%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equalize the channel impulse response by using a linear filter, sending
% known sequence and using LMS and RLS filters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
nsample = 10000; % sample size
refbit = randi([0 1], 2, nsample);
% mapping bit to symbol
sym = symbolizer_mqam(refbit);
% channel impulse response
h = [0.01, 0.2, 0.25, 0.5, 1.0] + 1i * [0.01, 0.2, 0.25, 0.5, 1.0];
% h = h ./ sum(h);

% x = conv(sym, h, 'same');
% x = fftfilt(h, sym);
% x = filter(h, 1, sym);

% filtering the input signal
x = zeros(size(sym));
p = length(h);
sym_ext = [zeros(p - 1, 1); sym];
for ii = 1 : nsample
    x(ii) = h * sym_ext((1 : p) + (ii - 1));
end
x = x(:);
sigma2 = calcrms(x)^2 / 50; % noise power
x = x + gaussian_noise(nsample, 1, sigma2, 'linear', 'complex');

yy_lms = least_squares_filter(x, sym, 'LMS', .01, [], 2*p);
err_lms = yy_lms - sym;

% in this example, smaller forgetting factor gives better result
yy_rls = least_squares_filter(x, sym, 'RLS', [], .55, 2*p);
err_rls = yy_rls - sym;

sym_dec = hard_decision(yy_rls, 4);
figure; plot(abs(sym - sym_dec), 'o'); grid on; ylim([-2 2]);

figure;
plot(real(err_lms)); hold on; 
plot(real(err_rls)); grid on;
ylim([-1, 1]);
xlabel('samples'); ylabel('error');
legend('LMS', 'RLS');

figure;
plot(yy_lms, 'b.'); hold on;
plot(yy_rls, 'r.'); grid on;
xlim([-2 2]); ylim([-2 2]);
