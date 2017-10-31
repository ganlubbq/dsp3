%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equalize the channel impulse response by using a linear filter, sending
% known sequence and using LMS and RLS filters.
%
% For static system, the forgetting factor of RLS could be 1.0. The smaller
% the forgetting factor, the shorter memory of the RLS algorithm and more
% capability of tracking time-varying system.
%
% Try to compare the filter taps convergence in cases of lambda = 1.0 and
% lambda = 0.5. In extreme case with very small lambda, RLS becomes almost
% memoryless, turnning "squares" to "square", i.e., focusing on the current
% sample only.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
nsample = 1000;

%%% complex signal
refbit = randi([0, 1], 2, nsample);
sym = symbolizer_mqam(refbit);

%%% try real signal
% sym = real(sym);

%%% channel impulse response
h = [-0.02, 0.5, 1.0];
% h = h ./ sum(h);

%%% filtering the input signal
% x = conv(sym, h, 'same');
% x = fftfilt(h, sym);
% x = filter(h, 1, sym);
x = zeros(size(sym));
p = length(h);
sym_ext = [zeros(p - 1, 1); sym];
for ii = 1 : nsample
    x(ii) = h * sym_ext((1 : p) + (ii - 1));
end
x = x(:);
sigma2 = calcrms(x)^2 / 100; % noise power
x = x + gaussian_noise(nsample, 1, sigma2, 'linear', 'complex');

[x_lms, w_lms] = least_squares_filter(x, sym, 'LMS', .01, [], 5);
[x_rls, w_rls] = least_squares_filter(x, sym, 'RLS', [], 1.0, 5);
err_lms = x_lms - sym;
err_rls = x_rls - sym;

%%% Compare bit error
% sym_dec = hard_decision(x_rls, 4);
% figure; plot(abs(sym - sym_dec), 'o'); grid on; ylim([-2 2]);

figure;
plot(real(err_lms)); hold on; 
plot(real(err_rls)); grid on;
ylim([-1, 1]);
hold off
xlabel('samples'); 
ylabel('error');
legend('LMS', 'RLS');

figure;
plot(x_lms, 'b.'); hold on;
plot(x_rls, 'r.'); grid on;
xlim([-3, 3]); 
ylim([-3, 3]);
hold off

figure;
plot(real(w_lms.'), 'Color', color_table(1), 'LineWidth', 6); hold on
plot(real(w_rls.'), 'Color', color_table(2), 'LineWidth', 2); grid on
hold off
xlabel('samples'); ylabel('filter taps');
