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
w = gaussian_noise(nsample, 1, sigma2, 'linear', 'real');
x = x + w;


% estimation initialization - lms
mu = 0.01;
he_lms = zeros(p, 1);
err_lms = zeros(nsample, 1);
for ii = 1 : nsample
    ss = s_ext((1 : p) + (ii - 1));
    he_lms = he_lms - mu * (transpose(he_lms) * ss - x(ii)) * conj(ss);
    err_lms(ii) = transpose(he_lms) * ss - x(ii);
end


% estimation initialization - rls
he_rls = zeros(p, 1);
Sigma = 1e5 * eye(p); % covariance matrix of estimation
err_rls = zeros(nsample, 1);
for ii = 1 : nsample
    ss = s_ext((1 : p) + (ii - 1));
    gain = Sigma * ss / (sigma2 + transpose(ss) * Sigma * ss);
    he_rls = he_rls + gain * (x(ii) - transpose(ss) * he_rls);
    Sigma = (eye(p) - gain * transpose(ss)) * Sigma;
    err_rls(ii) = x(ii) - transpose(ss) * he_rls;
end

figure;
plot(err_lms); hold on; 
plot(err_rls); grid on; 
ylim([-1, 1]);
xlabel('samples'); ylabel('error');
legend('LMS', 'RLS');
hold off;
