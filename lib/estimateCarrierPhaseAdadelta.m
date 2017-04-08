function [theta, s, J] = estimateCarrierPhaseAdadelta(signals, observations, theta_ini)
% Implementing the phase-locked loop to estimate the carrier phase based on
% the least squares criteria and adadelta algorithm (variable step size).
% Signals in the input is known as a reference.
if nargin < 3
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
gamma = 0.95;
beta  = 0.99;
meanSquareGradient = 0;
meanSquareDeltaPhi = 0;
epsilon = 1e-2;

for k = 2 : length(observations)
    % output symbols with phase correction
    s(k) = observations(k) .* exp(-1i * theta(k-1));
    
    % stochastic gradient
    grad(k) = -imag(s(k) .* conj(signals(k)));

    % update mean square of gradient using running sum with exponential
    % decay
    meanSquareGradient = gamma * meanSquareGradient + (1 - gamma) * grad(k)^2;
    
    % rms of gradient
    rmsGradient = sqrt(meanSquareGradient + epsilon);
    
    % rms of delta phi, kicked off by \epsilon
    rmsDeltaPhi = sqrt(meanSquareDeltaPhi + epsilon);
    
    % the adaptive delta
    deltaPhi(k) = (rmsDeltaPhi / rmsGradient) * grad(k);
    
    % update mean square of delta using running sum with exponential
    % deay
    meanSquareDeltaPhi = gamma * meanSquareDeltaPhi + (1 - gamma) * deltaPhi(k)^2;
    
    % delta integration using running sum with exponential decay
    gint(k) = beta * gint(k-1) + (1 - beta) * deltaPhi(k);
    
    % update filter coeff. along opposite direction of gradient
    theta(k) = theta(k-1) - deltaPhi(k) - gint(k);

    % squared error
    J(k) = abs(s(k) - signals(k)).^2;
end

return
