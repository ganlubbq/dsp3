% Get scale factor for scaling the canonical constellation to max{1,1} form
% 
% Example: getScaleFactorQAM( M )
% 
% Input: 
% 
% Reference: 
% 
% Note: Get a canonical form first and then devided by this factor would
% result in a max{1,1} form
% 
% See Also: getPowerFactorQAM
% 
% Copyright 2015 default

function [ s ] = getScaleFactorQAM( M )

s = sqrt(ALPHABET_SIZE)-1; % for square mQAM

if M==2
	s = 1;
elseif M==8
	s = 2;
elseif M==32
	s = 5;
else
	error('unsupported format');
end

return

