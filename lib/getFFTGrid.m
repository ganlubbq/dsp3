function freqGrid = getFFTGrid(nsample, fs)
% Get two-sided zero-first frequency grid vector for FFT
% 
% Example: freqGrid = getFFTGrid(nsample, fs)
% 
% Input: 
% 
% Note: 
% 
% See Also: fftshift

if mod(nsample, 2)
  freqGrid = [(0 : 1 : (nsample-1)/2), (-(nsample-1)/2 : 1 : -1)] * fs / nsample;
else
  freqGrid = [(0 : 1 : nsample/2-1), (-nsample/2 : 1 : -1)] * fs / nsample;
end

return
