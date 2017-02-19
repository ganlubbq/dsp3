function [theta, s, J] = estimateCarrierPhaseLMS_hardDecision(observations, mn, stepsize, theta_ini)
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
else
	mu1 = stepsize(1);
	mu2 = stepsize(2);
end


for k = 2:length(observations)
    
    % output
    s(k) = observations(k) .* exp(-1j * theta(k-1));
    
    signal = slicer(s, mn);
    
    % stochastic gradient
    grad(k) = -imag(s(k) .* conj(signal));
    
    % error integration
    eint(k) = eint(k-1) + grad(k);
    
    % update filter coeff. along opposite direction of gradient
    theta(k) = theta(k-1) - mu1*grad(k) - mu2*eint(k);
    
    % squared error
    J(k) = abs(s(k) - signals(k)).^2;
end

return

% slice columns of complex input into symbols 
function d = slicer(x,mn)

xnorm = normalizeQam(x(:),mn);

switch mn % MSB first
    case 2
        bit = msign(real(xnorm));
        
	case 4
		bit = [msign(real(xnorm)), ~msign(imag(xnorm))];
		
	case 16
		bound = 2;
		bit = [msign(real(xnorm)), msign(real(xnorm)-ksign(real(xnorm))*bound),...
			~msign(imag(xnorm)), ~msign(imag(xnorm)-ksign(imag(xnorm))*bound)];
		
	case 64 % 64qam slicing is combination of qpsk and 16qam, i.e.
		% [4qam 16qam 16qam 4qam 16qam 16qam]
		bound = 2;
		
		% 4qam header
		bit1 = [msign(real(xnorm)), ~msign(imag(xnorm))];
		
		% convert 64qam symbol in one quater to 16qam
		x = x-(ksign(real(xnorm))*4+ksign(imag(xnorm))*4j);
		
		% slicing 16qam
		bit2 = [msign(real(xnorm)), msign(real(xnorm)-ksign(real(xnorm))*bound),...
			~msign(imag(xnorm)), ~msign(imag(xnorm)-ksign(imag(xnorm))*bound)];
		
		bit = [bit1(:,1) bit2(:,1:2) bit1(:,2) bit2(:,3:4)];
		
	case 256
		% TBA
        
	otherwise
		error('unsupported format');
end

% convert bit in rows to dec by MSB to LSB order
twos = 2.^(log2(mn)-1:-1:0);
twos = ones(length(xnorm),1) * twos;

% index of symbols from topleft to bottomright by columns
ndx = sum(bit.*twos,2) + 1;

% symbols from topleft to bottomright by columns
c = constellation(mn);

d = c(ndx);

return

% return 1 if x>0 and 0 if x<=0
function s = msign(x)

s = (x>0);

return

% return 1 if x>0 and -1 if x<=0
function s = ksign(x)

s = (x>0)*2-1;

return

