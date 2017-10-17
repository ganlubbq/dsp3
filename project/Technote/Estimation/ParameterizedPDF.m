% This script plots the parameterized pdf of classical estimation model,
% where the unknown parameter is considered to be deterministic. Note that
% slicing along the fixed data would generate a likelihood function of
% unknown parameter.
%
% Kay, Steven M. "Fundamentals of statistical signal processing: estimation
% theory." (1993).
%
clear

% nominal data
x = -1 : 0.05 : 1;

% the deterministic unknown parameter, gaussian mean
theta = -1 : 0.05 : 1;

% PDF = 1/sqrt(2*pi*0.1)*exp(-(x.-theta).^2/(2*0.1));
% PDF = PDF';

% varying with data
PDF3 = [];
for jj = 1 : length(theta)
    PDF = 1/sqrt(2*pi*0.1) * exp(-(x-theta(jj)).^2/(2*0.1));
    PDF = PDF(:);
    PDF3 = [PDF3, PDF];
end

% plot the parameterized pdf
figure;
subplot(211); mesh(theta,x,PDF3); view([-64 19]);
xlabel('x'); ylabel('\theta'); zlabel('PDF');

% the likelihood function can be obtained by fix data x and varying unknown
% parameter
subplot(212);
plot(theta, PDF3(10,:)); 
xlabel('\theta'); ylabel('p(\theta|x)');
grid on; box on
