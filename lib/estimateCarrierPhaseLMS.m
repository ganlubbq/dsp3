function [theta, s, J] = estimateCarrierPhaseLMS(signals, observations, stepsize, theta_ini)
%ESTIMATECARRIERPHASELMS Implementing the adaptive filter to estimate the
%carrier phase based on the Least Squares criteria. Signals in the input is
%known as a reference. Schocastic gradient descent is applied. Loop filter
%up to 2nd-order is implemented

if nargin < 4
    theta_ini = 0;
end

s = zeros(size(observations));
J = zeros(size(observations));
grad = zeros(size(observations));
eint = zeros(size(observations));
theta = zeros(size(observations));

% initialize stochastic gradient algorithm with least squares criteria
s(1) = observations(1);
theta(1) = theta_ini;
eint(1) = 0;


if length(stepsize) == 1
	mu1 = stepsize;
    mu2 = 0;
else
	mu1 = stepsize(1);
	mu2 = stepsize(2);
end


for k = 2:length(observations)
    
    % output
    s(k) = observations(k) .* exp(-1j * theta(k-1));
    
    % stochastic gradient
    grad(k) = -imag(s(k) .* conj(signals(k)));
    
    % error integration
    eint(k) = eint(k-1) + grad(k);
    
    % update filter coeff. along opposite direction of gradient
    theta(k) = theta(k-1) - mu1*grad(k) - mu2*eint(k);
    
    % squared error
    J(k) = abs(s(k) - signals(k)).^2;
end

return
