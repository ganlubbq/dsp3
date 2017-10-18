%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chracterize (estimate) the channel impulse response by sending known 
% sequence and using LMS and RLS filters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

nsample = 2000; % sample size

s = 2 * (randi([0 1], nsample, 1) - 0.5);

h = [0.01, 0.2, 0.25, 0.5, 1.0];
h = h ./ sum(h);

%c1 = conv(s, h, 'same');
%c2 = fftfilt(h, s);
c3 = filter(h, 1, s);

x = zeros(size(c3));
p = length(h);
for ii = 1 : nsample - p
    ss = s((1:p) + (ii-1));
    x(ii) = fliplr(h) * ss;
end
x = x(:) * 10;

sigma2 = .1; % noise power
w = gaussian_noise(nsample, 1, sigma2, 'linear', 'real');

x = x + w;


% estimation initialization - lms
mu = 0.01;
he_lms = zeros(p, 1);

for ii = 1 : nsample - p
    ss = s((1:p) + (ii-1));
    he_lms = he_lms - mu * (transpose(he_lms) * ss - x(ii)) * conj(ss);
    err_lms(ii) = transpose(he_lms) * ss - x(ii);
end

% estimation initialization - rls
he_rls = zeros(p, 1);
Sigma = 100 * eye(p); % covariance matrix of estimation

for ii = 1 : nsample - p
    ss = s((1:p) + (ii-1));
    gain = Sigma * ss / (1 + transpose(ss) * Sigma * ss);
    he_rls = he_rls + gain * (x(ii) - transpose(ss) * he_rls);
    err_rls(ii) = x(ii) - transpose(ss) * he_rls;
end

figure;
subplot(211); plot(err_lms); grid on; xlabel('samples'); ylabel('error');
legend('LMS');
subplot(212); plot(err_rls); grid on; xlabel('samples'); ylabel('error');
legend('RLS');
