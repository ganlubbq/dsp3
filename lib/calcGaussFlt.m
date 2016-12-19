% Calculate frequency response of electrical super Gaussian filter
% 
% Example: [H] = calcGaussFlt(nSample,fs,fo,order,bandwidth)
% 
% Input: 
%       nSample     - number of samples
%       fs          - sampling frequency [Hz]
%       fo          - offset frequency [Hz]
%       order       - order of filter
%       bandwidth   - 3dB bandwidth [Hz] one-sided
% 
% Reference: Martin Pfennigbauer et al "Choice of MUX/DEMUX filter
% characteristics for NRZ RZ and CSRZ..." JLT 2006
% 
% Note: 
% 
% See Also: calcOptGaussFlt
% 
% Copyright 2015 Dawei Wang 
% 
% Email: dawei.zju@gmail.com 

function [H] = calcGaussFlt(nSample,fs,fo,order,bandwidth)

% frequency interval of samples
deltaFs = fs/nSample;

% frequency grid in Hz
freqGrid = [(0:nSample/2-1)'; flipud(-(1:(nSample/2))')] *deltaFs;

% frequency response, amplitude spetrum
H = exp(-0.5*log(2)*((freqGrid-fo)/bandwidth).^(2*order));

return

