function theta = estimateCarrierPhaseLMS(signals,observations,stepsize)
%ESTIMATECARRIERPHASELMS Implementing the adaptive filter to estimate the
%carrier phase based on the Least Squares criteria. Signals in the input is
%known as a reference. Schocastic gradient descent is applied.

% initialize stochastic gradient algorithm with least squares criteria
s(1) = observations(1);
theta(1) = 0;

for k = 2:length(observations)
    
    % output
    s(k) = observations(k).*exp(-1j*theta(k-1));
    
    % stochastic gradient, which is equivalent to sin(theta), a typical
    % S-Curve for phase error detector
    grad(k) = -imag(s(k).*conj(signals(k)));
    
    % update filter coeff. along opposite direction of gradient
    theta(k) = theta(k-1) - stepsize*grad(k);
    
    % squared error
    J(k) = abs(s(k)-signals(k)).^2;
end

return


