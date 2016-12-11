% Get frequency grid vector for FFT
% 
% Example: freqGrid = getFFTGrid(nSample,fs)
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

function freqGrid = getFFTGrid(nSample,fs)

% frequency grid in Hz
freqGrid = [(0:nSample/2-1)'; flipud(-(1:(nSample/2))')] *fs/nSample;

return

