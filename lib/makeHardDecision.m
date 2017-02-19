function [symbol] = makeHardDecision(x, mn)
%MAKEDECISION % slice columns of input into symbols 
%
%   Copyright 2012

if ~iscolumn(x)
    error('first input has to be a column vector');
end

xnorm = normalizeQam(x, mn);

switch mn % MSB first
    case 2
        bit = msign(real(xnorm));
        
	case 4
		bit = [msign(real(xnorm)), ~msign(imag(xnorm))];
		
	case 16
		bound = 2;
		bit = [msign(real(xnorm)), msign(real(xnorm)-ksign(real(xnorm))*bound),...
			~msign(imag(xnorm)), ~msign(imag(xnorm)-ksign(imag(xnorm))*bound)];
		
	case 64 
        % 64qam slicing is combination of qpsk and 16qam, i.e.
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
		error('unsupported modulation format');
end

% convert bit in rows to dec by MSB to LSB order
twos = 2.^(log2(mn)-1:-1:0);
twos = ones(length(xnorm),1) * twos;

% index of symbols from topleft to bottomright by columns
ndx = sum(bit.*twos,2) + 1;

% symbols from topleft to bottomright by columns
c = constellation(mn);

symbol = c(ndx);

return

% return 1 if x>0 and 0 if x<=0
function s = msign(x)

s = (x>0);

return

% return 1 if x>0 and -1 if x<=0
function s = ksign(x)

s = (x>0)*2-1;

return

