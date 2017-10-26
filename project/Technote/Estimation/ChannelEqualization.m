%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equalize the channel impulse response by using a linear filter, sending
% known sequence and using LMS and RLS filters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

nsample = 4000; % sample size
s = 2 * (randi([0, 1], nsample, 1) - 0.5);
h = [0.01, 0.2, 0.25, 0.5, 1.0];
h = h ./ sum(h);
p = length(h);

%x = conv(s, h, 'same');
%x = fftfilt(h, s);
%x = filter(h, 1, s);

x = zeros(size(s));
s_ext = [zeros(p - 1, 1); s];
for ii = 1 : nsample
    ss = s_ext((1 : p) + (ii - 1));
    x(ii) = h * ss;
end
x = x(:);
sigma2 = calcrms(x)^2 / 100; % noise power
x = x + gaussian_noise(nsample, 1, sigma2, 'linear', 'real');

[x_lms, w_lms] = least_squares_filter(x, s, 'LMS', .01, [], 2*p);
err_lms = x_lms - s;
[x_rls, w_rls] = least_squares_filter(x, s, 'RLS', [], 1.0, 2*p);
err_rls = x_rls - s;

figure;
plot(err_lms); hold on; 
plot(err_rls); grid on;
hold off
ylim([-1, 1]);
xlabel('samples'); ylabel('error');
legend('LMS', 'RLS');

figure;
plot(w_lms.'); hold on
plot(w_rls.'); grid on
hold off
