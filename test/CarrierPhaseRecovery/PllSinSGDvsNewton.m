% TEST SCRIPT FOR IMPLEMENTING A 2ND-ORDER PHASE-LOCK LOOP BASED ON SGD
% ALGORITHM
%
% SGD means stochastic gradient descent
%
% Note that multiple samples per symbol are used for the loop. In the case
% of large carrier frequency offset, it is necessary to turn on the
% 2nd-order loop filter. One can turn on and off the 2nd order loop filter
% to observe the phase tracking results.
%
% Kay, Steven M. "Fundamentals of statistical signal processing: estimation
% theory." (1993).
%
% SINGLE FREQUENCY WITH FIXED PHASE ERROR
function [] = PllSinSGDvsNewton()

frequency = 10;

% sampling frequency
fs = 80;

tvec = 0 : (1/fs) : 20;

% single frequency with phase error
x = exp(1i * (2*pi*frequency*tvec + pi*2.01/4));
n = genWGN(size(x,1),size(x,2),0.002,'linear','complex');
x = x + n;

% reference frequency with offset
cfo = frequency * 0.06;
ref = exp(1i * (2*pi*(frequency + cfo)*tvec));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the nonlinear LS function with gradient descent method with data
% model J = |x*exp(-i*phi) - ref|^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% _Using multiple samples per symbol_
mu1 = 0.05;         % gain parameter 1st order
mu2 = 0.001 * 1;    % gain parameter 2nd order
s(1) = x(1); 
phi(1) = 0; 
nco(1) = 0;
for k = 2:length(x)
    % output
    s(k) = x(k) .* exp(-1i * phi(k-1));
    
    % stochastic gradient; PED with s-curve, equivalent to sin(angle())
    grad(k) = - imag(s(k) .* conj(ref(k)));
    
    % err integration
    nco(k) = nco(k-1) + grad(k);
    
    % update filter coeff. along opposite direction of gradient, equivalent
    % to a low-pass filter
    phi(k) = phi(k-1) - mu1*grad(k) - mu2*nco(k);
    
    % squared error
    J(k) = abs(s(k) - ref(k)).^2;
end

h1 = figure; title('gradient decent');
% phase estimation
subplot(211); plot(tvec, mod(phi, 2*pi)); grid on
% Squares
subplot(212); plot(dbw(J(1 : 1000))); grid on; ylim([-100 20])
title('SGD');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving transformed linear least squares function with gradient descent
% method with data model J = |ref*theta - x|^2, by viewing \theta =
% exp(i*phi) as the unknown parameter. However, this simple transformation
% results a NONLINEAR constrained LS problem...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% no solution yet...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the nonlinear least squares function with Newton-Raphson method
% the convergence problem can be observed when the system has larger noise
s = [];
phi = [];
grad = [];
s(1) = x(1);
phi(1) = 0;

% add this to hessian can solve the convergence problem to some extent
epsilon = 0.5; 

for k = 2:length(x)
    % output
    s(k) = x(k) .* exp(-1i * phi(k-1));
    
    % stochastic gradient
    grad(k) = - imag(s(k) .* conj(ref(k)));
    
    % stochastic hessian (BE CAREFUL when hessian is small!!!)
    hess(k) = real(s(k) .* conj(ref(k))) + epsilon;
    
    % solving the next root of linearization model
    phi(k) = phi(k-1) - grad(k) / hess(k);
    
    % squared error
    J(k) = abs(s(k) - ref(k)).^2;
end

h2 = figure; 
% phase estimation
subplot(211); plot(tvec, mod(phi, 2*pi)); grid on
% Squares
subplot(212); plot(dbw(J(1 : 1000))); grid on; ylim([-100 20])
title('Newton-Raphson');

% depending on the location of initial guess, this method may converge to
% another solution with pi shift relatively

mngFigureWindow(h1,h2);

return
