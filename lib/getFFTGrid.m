function freq = getFFTGrid(nsample, fs)
% Get two-sided zero-first frequency grid vector for FFT
% 
% Example: freq = getFFTGrid(nsample, fs)
% 
% Input: 
% 
% Note: For even and odd number of samples, the frequency vector centers at
% zero. Even number of samples gets one extra point at the negative end.
% 
% See Also: fftshift

if mod(nsample, 2)
  freq = [(0 : 1 : (nsample - 1)/2), (-(nsample - 1)/2 : 1 : -1)] * fs / nsample;
else
  freq = [(0 : 1 : nsample/2 - 1), (-nsample/2 : 1 : -1)] * fs / nsample;
end

return
