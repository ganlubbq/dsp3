function [theta] = estimateCarrierPhaseVV(x,block_size,avg_mode)
% Feedforward carrier phase estimation using VVPE.
%
%			avg_mode	- 0 for block 1 for sliding window
%   
% See also:

% raise the phase to 4th power not the amplitude
x4 = abs(x).*exp( 1j*angle( x.^ 4 ) );

if avg_mode==0
	% make sure length of x4 can be divided by block size
	tmp = mod(length(x4),block_size);
	if tmp
		x4 = [x4(:);zeros(tmp,1)];
	end
	x4b = reshape(x4(:),block_size,[]);
	x4s = sum(x4b)/block_size;	
elseif avg_mode==1
	% sliding average
	x4s = smooth(x4(:),block_size);
else
	error('incorrect mode number');
end

% unwrapping
theta = unwrap(angle(x4s))/4;

if avg_mode==0
	theta = theta(:)*ones(1,block_size);
	theta = theta.';
	theta = theta(1:end);
end

return


