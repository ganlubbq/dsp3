function photocurrent = photoDetectorSE(Ex, responsivity)
% DESCRIPTION
% 
% Example: photocurrent = photoDetectorSE(Ex, responsivity)
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: photoDetectorBN

photocurrent = responsivity .* abs(Ex).^2;

return