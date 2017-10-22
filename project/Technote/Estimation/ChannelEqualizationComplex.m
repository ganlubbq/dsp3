%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equalize the channel impulse response by using a linear filter, sending
% known sequence and using LMS and RLS filters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
nsample = 12000; % sample size
refbit = randi([0 1], 2, nsample);
% mapping bit to symbol
sym = symbolizer_mqam(refbit);
% channel impulse response
h = [0.01, 0.2, 0.25, 0.5, 1.0];
h = h ./ sum(h);

%c1 = conv(sym, h, 'same');
%c2 = fftfilt(h, sym);
c3 = filter(h, 1, sym);

% filtering the input signal
x = zeros(size(c3));
p = length(h);
for ii = 1 : nsample - p
    ss_real = real(sym((1 : p) + (ii - 1)));
    ss_imag = imag(sym((1 : p) + (ii - 1)));
    x(ii) = fliplr(h) * ss_real + 1i * fliplr(h) * ss_imag;
end
x = x(:);
sigma2 = .05; % noise power
w = gaussian_noise(nsample, 1, sigma2, 'linear', 'complex');
x = x + w;

% linear estimation initialization - lms
mu = 0.01;
he_lms = zeros(p, 1); he_lms(3) = 1;
err_lms = zeros(p, 1);
yy_lms = zeros(p, 1);
for ii = 1 : nsample - p
    xx = x((1 : p) + (ii - 1));
    he_lms = he_lms - mu * (he_lms.' * xx - sym(ii)) * conj(xx);
    err_lms(ii) = he_lms.' * xx - sym(ii);
    yy_lms(ii) = he_lms.' * xx;
end

% estimation initialization - rls
he_rls = zeros(p, 1);
Sigma = 10 * eye(p); % covariance matrix of estimation
err_rls = zeros(p, 1);
yy_rls = zeros(p, 1);
for ii = 1 : nsample - p
    xx = x((1 : p) + (ii - 1));
    gain = Sigma * xx / (1 + xx' * Sigma * xx);
    he_rls = he_rls + gain * (sym(ii) - xx' * he_rls);
    err_rls(ii) = sym(ii) - xx' * he_rls;
    yy_rls(ii) = he_rls' * xx;
end

figure;
plot(abs(err_lms)); hold on; xlabel('samples'); ylabel('error');
plot(abs(err_rls)); grid on;
legend('LMS', 'RLS');

figure;
plot(x, '.'); hold on;
plot(yy_lms, 'x');
plot(yy_rls, 'o'); grid on;
xlim([-2 2]);
ylim([-2 2]);
