function [H] = calcOptGaussFlt(nSample,fs,fo,order,bandwidth)
% Calculate frequency response of optical super Gaussian filter
% 
% Example: [H] = calcOptGaussFlt(nSample,fs,fo,order,bandwidth)
% 
% Input: 
%       nSample     - number of samples
%       fs          - sampling frequency [Hz]
%       fo          - offset frequency [Hz]
%       order       - order of filter
%       bandwidth   - 3dB bandwidth [Hz] two-sided
% 
% Reference: Martin Pfennigbauer et al "Choice of MUX/DEMUX filter
% characteristics for NRZ RZ and CSRZ..." JLT 2006
% 
% Note: 
% 
% See Also: 
% 
% Copyright 2015 Default

% frequency interval of samples
deltaFs = fs/nSample;

% frequency grid in Hz
freqGrid = [(0:nSample/2-1)'; flipud(-(1:(nSample/2))')] *deltaFs;

% frequency response, amplitude spetrum
H = exp(-0.5*log(2)*(2*(freqGrid-fo)/bandwidth).^(2*order));

return
