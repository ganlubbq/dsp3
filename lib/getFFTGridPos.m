function freqGrid = getFFTGridPos(nsample, fs)
% Get one-sided positive only frequency grid vector for FFT
% 
% Example: freqGrid = getFFTGridPos(nsample, fs)
% 
% Input: 
% 
% Note: 
% 
% See Also: fftshift

if mod(nsample, 2)
  freqGrid = (0 : 1 : (nsample-1)/2) * fs / nsample;
else
  freqGrid = (0 : 1 : nsample/2-1) * fs / nsample;
end

return
