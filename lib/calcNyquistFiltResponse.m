function H = calcNyquistFiltResponse(nSample,fs,fbaud,alpha,mode)
% Calculate frequency response of Nyquist filter with input roll-off and
% baud rate
% 
% Example: 
%   H = calcNyquistFiltResponse(nSample,fs,fbaud,alpha,mode)
% 
% Input: 
%       nSample     - number of samples
%       fs          - sampling frequency [Hz]
%       order       - filter order
%       bandwidth   - filter 3dB bandwidth in [Hz]
% 
% Reference: 
% 
% Note: this filter calculate a raised cosine filter first, then combined
% with a reversed sinc filter
% 
% See also: calcDispResponse, calcOptGaussFlt, calcRcosResponse
% 
% Copyright 2015 Default

HRCOS = calcRCFreqResponse(nSample,fs,fbaud,alpha,mode);

freq = [(0:nSample/2-1)'; flipud(-(1:(nSample/2))')] * fs / nSample;

f_high = (1+alpha)*fbaud/2;

HSINC = ones(size(HRCOS));
HSINC(abs(freq) <= f_high) = sinc(freq(abs(freq) <= f_high) ./ fbaud);

H = HRCOS./HSINC;

return
