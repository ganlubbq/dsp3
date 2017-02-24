function freqGrid = getFFTGridPos(nSample, fs)

% frequency interval of samples
deltaFs = fs / nSample;

% frequency grid in Hz
freqGrid = (0:nSample/2-1)' * deltaFs;

return