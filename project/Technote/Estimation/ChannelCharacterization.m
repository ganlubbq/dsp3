%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chracterize (estimate) the channel impulse response by sending known 
% sequence and using LMS and RLS filters.
%
% For stationary process, the forgetting factor of RLS could be 1.0. The
% smaller the forgetting factor, the shorter memory of the RLS algorithm. 
%
% Try to compare the filter taps convergence in cases of lambda = 1.0 and
% lambda = 0.5. In extreme case with very small lambda, RLS becomes almost
% memoryless, turnning "squares" to "square", i.e., focusing on the current
% sample only.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

nsample = 2000; % sample size
s = 2 * (randi([0 1], nsample, 1) - 0.5);
h = [0.01, 0.2, 0.25, 0.5, 1.0];
h = h(:) ./ sum(h);
p = length(h);

%x = conv(s, h, 'same');
%x = fftfilt(h, s);
%x = filter(h, 1, s);

x = zeros(size(s));
s_ext = zero_pad(s, p-1, 'front');
for ii = 1 : nsample
    ss = s_ext((1 : p) + (ii - 1));
    x(ii) = h.' * ss;
end
x = x(:);
sigma2 = calcrms(x)^2 / 100; % noise power
x = x + gaussian_noise(nsample, 1, sigma2, 'linear', 'real');

[y_lms, h_lms] = least_squares_filter(s, x, 'BPSK', 'LMS', .01, [], p);
[y_rls, h_rls] = least_squares_filter(s, x, 'BPSK', 'RLS', [], 1.0, p);
err_lms = x - y_lms;
err_rls = x - y_rls;

figure;
plot(err_lms); hold on; 
plot(err_rls); grid on; 
ylim([-1, 1]);
xlabel('samples'); ylabel('error');
legend('LMS', 'RLS');
hold off;

figure;
plot(transpose(h * ones(1, size(h_lms, 2))), 'k'); hold on
plot(h_lms.', 'Color', color_table(1), 'LineWidth', 4); grid on
plot(h_rls.', 'Color', color_table(2), 'LineWidth', 2); grid on
xlabel('samples'); ylabel('filter taps');
ylim([-0.1, 0.6]);
hold off

