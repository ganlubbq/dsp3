% TEST SCRIPT FOR IMPLEMENT A 2ND-ORDER PHASE-LOCK LOOP BASED ON SGD ALGORITHM
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
function [] = PllSineSGDvsNewton()

frequency = 10;

% sampling frequency
fs = 80;

tvec = 0:(1/fs):20;

% single frequency with phase error
x = exp(1j*(2*pi*frequency*tvec+pi*2.01/4));
n = genWGN(size(x,1),size(x,2),0.001,'linear','complex');
x = x + n;

% reference frequency with offset
cfo = 0.5;
ref = exp(1j*(2*pi*(frequency+cfo)*tvec));

%% Solving the nonlinear LS function with gradient descent method
% Using least squares equalization model, i.e., J = |x*exp(-j*phi)-ref|^2

% initialize stochastic gradient descent algorithm, implemented as a
% phase-lock loop. 
% _Using multiple samples per symbol_
mu1 = 0.05;  % gain parameter 1st order
mu2 = 0.001 * 1; % gain parameter 2nd order
s(1) = x(1);  % output
phi(1) = 0;  % phase estimation
nco(1) = 0;
for k = 2:length(x)
    % output
    s(k) = x(k) .* exp(-1i * phi(k-1));
    
    % stochastic gradient; PED with s-curve, equivalent to sin(angle())
    grad(k) = -imag(s(k) .* conj(ref(k)));
    
    % err integration
    nco(k) = nco(k-1) + grad(k);
    
    % update filter coeff. along opposite direction of gradient, equivalent
    % to a low-pass filter
    phi(k) = phi(k-1) - mu1*grad(k) - mu2*nco(k);
    
    % squared error
    J(k) = abs(s(k) - ref(k)).^2;
end

h1 = figure; title('gradient decent');
% show the phase difference
subplot(2,2,1); plot(tvec,real(x),'b',tvec,real(ref),'r'); grid on
% compare reference with output
subplot(2,2,2); plot(tvec,real(ref),'b',tvec,real(s),'r'); grid on
% phase estimation
subplot(2,2,3); plot(tvec,phi); grid on
% Squares
subplot(2,2,4); plot(tvec,dbw(J)); grid on; ylim([-100 20])

% Solving transformed linear LS function with gradient descent method
% Using least squares data model, i.e., J = |ref*theta-x|^2, i.e., viewing
% theta=exp(i*phi) as the unknown parameter. However, this simple
% transformation results a NONLINEAR constrained LS problem...

% no solution yet...

% Solving the nonlinear LS function with Newton-Raphson method
% Using least squares equalization model, i.e., J = |x*exp(-j*phi)-ref|^2
s = [];
phi = [];
grad = [];
% initialize stochastic Newton-Raphson algorithm, the NR algorithm suffers
% from convergence problem for particular initial points.
mu1 = 0.05;  % gain parameter 1st order
mu2 = 0.001 * 1; % gain parameter 2nd order
s(1) = x(1);
phi(1) = 0;
nco(1) = 0;
for k = 2:length(x)
    % output
    s(k) = x(k) .* exp(-1i * phi(k-1));
    
    % stochastic gradient
    grad(k) = -imag(s(k) .* conj(ref(k)));
    
    % stochastic hessian (BE CAREFUL when hessian is small!!!)
    hess(k) = real(s(k) .* conj(ref(k)));
    
    % err integration
    nco(k) = nco(k-1) + grad(k) / hess(k);
    
    % solving the next root of linearization model
    phi(k) = phi(k-1) - mu1*grad(k)/hess(k) - mu2*nco(k);
    
    % squared error
    J(k) = abs(s(k) - ref(k)).^2;
end

h2 = figure; title('Newton-Raphson');
% show the phase difference
subplot(2,2,1); plot(tvec,real(x),'b',tvec,real(ref),'r'); grid on
% compare reference with output
subplot(2,2,2); plot(tvec,real(ref),'b',tvec,real(s),'r'); grid on
% phase estimation
subplot(2,2,3); plot(tvec,phi); grid on
% Squares
subplot(2,2,4); plot(tvec,dbw(J)); grid on; ylim([-100 20])

% depending on the location of initial guess, this method may converge to
% another solution with pi shift relatively

mngFigureWindow(h1,h2);

return
