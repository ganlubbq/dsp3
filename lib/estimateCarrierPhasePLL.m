function [theta, s, J] = estimateCarrierPhasePLL(signals, observations, stepsize, theta_ini)
% Implementing the phase-locked loop to estimate the carrier phase based on
% the least squares criteria and stochastic gradient descent algorithm.
% Signals in the input is known as a reference.
if nargin < 4
    theta_ini = 0;
end

s     = zeros(size(observations));
J     = zeros(size(observations));
grad  = zeros(size(observations));
gint  = zeros(size(observations));
theta = zeros(size(observations));

% initialize stochastic gradient algorithm with least squares criteria
s(1) = observations(1);
theta(1) = theta_ini;
gint(1) = 0;

if numel(stepsize) == 1
	mu1 = stepsize;
    mu2 = 0;
else
	mu1 = stepsize(1);
	mu2 = stepsize(2);
end

for k = 2 : length(observations)
    % output symbols with phase correction
    s(k) = observations(k) .* exp(-1i * theta(k-1));
    
    % stochastic gradient
    grad(k) = -imag(s(k) .* conj(signals(k)));
    
    % gradient integration
    gint(k) = gint(k-1) + grad(k);
    
    % update filter coeff. along opposite direction of gradient
    theta(k) = theta(k-1) - mu1 * grad(k) - mu2 * gint(k);
    
    % squared error
    J(k) = abs(s(k) - signals(k)).^2;
end

return
