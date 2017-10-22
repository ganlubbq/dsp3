%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the mean of WGN using LMS filter and RLS filter. The two filters
% can be considered as stochastic version of gradien descent and Gauss
% method, respectively.
% 
% The mean of WGN is assumed to be deterministic and hence this is classic
% estimation.
% 
% As expected, RLS filter converges much faster than LMS filter at expense
% of high computation complexity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

A = 10; % the true mean
nsample = 1000; % sample size
sigma2 = 1; % noise power
w = gaussian_noise(nsample, 1, sigma2, 'linear', 'real');
x = A + w;


% estimation initialization
mu = 0.01;
xe_lms = zeros(size(x));

gain = zeros(size(x)); % gain
xe_rls = zeros(size(x));
variance = zeros(size(x)); % variance of estimation

for ii = 1 : nsample
    if ii == 1
        xe_lms(ii) = x(1);
        xe_rls(ii) = x(1);
        variance(ii) = sigma2;
    else
        % LMS
        xe_lms(ii) = xe_lms(ii-1) + mu * (x(ii) - xe_lms(ii-1));
        % RLS
        gain(ii) = variance(ii-1) / (variance(ii-1) + sigma2);
        xe_rls(ii) = xe_rls(ii-1) + gain(ii) * (x(ii) - xe_rls(ii-1));
        variance(ii) = (1 - gain(ii)) * variance(ii-1);
    end
end

figure;
subplot(311); plot(dbw(variance)); grid on; xlabel('samples'); ylabel('variance');
subplot(312); plot(dbw(gain)); grid on; xlabel('samples'); ylabel('gain');
subplot(313); plot(dbw(xe_rls)); grid on; hold on; plot(xe_lms); 
xlabel('samples'); ylabel('estimated mean'); legend('RLS', 'LMS');
