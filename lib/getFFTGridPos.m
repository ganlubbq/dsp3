% DESCRIPTION
% 
% Example: 
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright 2015 Dawei Wang 
% 
% Email: dawei.zju@gmail.com 

function freqGrid = getFFTGridPos(nSample,fs)

% frequency interval of samples
deltaFs = fs/nSample;

% frequency grid in Hz
freqGrid = (0:nSample/2-1)' *deltaFs;

return

