function y = normalizeQam(x, mn)
% Normalize m-QAM signal to its canonical form
% 
% Example: y = normalizeQam (x, mn)
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: 

% get the scaling factor for canonical form, i.e., the average power of canonical
% form. Note that the scaling factor is 2/3*(mn-1) for square QAM
scale_factor = getPowerFactorQAM(mn);

% first normalize signal to UNIT average symbol energy
% then, multiply by scaling factor
y = x / sqrt(mean(abs(x).^2)) * sqrt(scale_factor);

return

function [p] = getPowerFactorQAM(M)

if M == 2
	p = 1;
elseif M == 8
	p = 10;
elseif M == 32
	p = 20;
else % for square mQAM
	p = 2/3 * (M - 1); 
end

return
