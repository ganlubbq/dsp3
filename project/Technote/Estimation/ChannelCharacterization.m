%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chracterize (estimate) the channel impulse response by sending known 
% sequence and using LMS and RLS filters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

nsample = 2000; % sample size
s = 2 * (randi([0 1], nsample, 1) - 0.5);
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

err_lms = x - least_squares_filter(s, x, 'LMS', .01, [], p);
err_rls = x - least_squares_filter(s, x, 'RLS', [], .9, p);

figure;
plot(err_lms); hold on; 
plot(err_rls); grid on; 
ylim([-1, 1]);
xlabel('samples'); ylabel('error');
legend('LMS', 'RLS');
hold off;
