function photocurrent = photoDetectorBN(Ep, En, Rp, Rn)
% DESCRIPTION
% 
% Example: photocurrent = photoDetectorBN(Ep, En, Rp, Rn)
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: photoDetectorSE

photocurrent = Rp .* abs(Ep).^2 - Rn .* abs(En).^2;

return
