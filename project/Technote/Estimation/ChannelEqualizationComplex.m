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
h = [0.01, 0.2, 0.25, 0.5, 1.0]% + 1i * [0.03, 0.05, 0.025, 0.02, 0.01];
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
sigma2 = calcrms(x)^2 / 100; % noise power
w = gaussian_noise(nsample, 1, sigma2, 'linear', 'complex');
x = x + w;

% linear estimation initialization - lms
mu = .01;
he_lms = zeros(p, 1);
err_lms = zeros(p, 1);
yy_lms = zeros(p, 1);
x_ext = [zeros(p - 1, 1); x];
for ii = 1 : nsample
    xx = x_ext((1 : p) + (ii - 1));
    he_lms = he_lms - mu * (xx.' * he_lms - sym(ii)) * conj(xx);
    err_lms(ii) = he_lms.' * xx - sym(ii);
    yy_lms(ii) = he_lms.' * xx;
end

% estimation initialization - rls
he_rls = zeros(p, 1);
Sigma = 1e5 * eye(p); % covariance matrix of estimation
err_rls = zeros(p, 1);
yy_rls = zeros(p, 1);
x_ext = [zeros(p - 1, 1); x];
for ii = 1 : nsample
    xx = x_ext((1 : p) + (ii - 1));
    gain = Sigma * xx / (sigma2 + xx.' * Sigma * xx);
    he_rls = he_rls + gain * (sym(ii) - xx.' * he_rls);
    Sigma = (eye(p) - gain * xx.') * Sigma;
    err_rls(ii) = sym(ii) - xx.' * he_rls;
    yy_rls(ii) = he_rls.' * xx;
end

figure;
plot(abs(err_lms)); hold on; 
plot(abs(err_rls)); grid on;
xlabel('samples'); ylabel('error');
legend('LMS', 'RLS');

figure;
plot(yy_lms, 'b.'); hold on;
plot(yy_rls, 'r.'); grid on;
xlim([-2 2]);
ylim([-2 2]);
