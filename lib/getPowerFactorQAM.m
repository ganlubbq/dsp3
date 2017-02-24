function [ p ] = getPowerFactorQAM( M )
% Get the power of canonical constellations for scaling
% 
% Example: 
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: normalize_mqam
% 
% Copyright 2015 

if M == 2
	p = 1;
elseif M == 8
	p = 10;
elseif M == 32
	p = 20;
else
	p = 2/3*(M-1); % for square mQAM
end

return