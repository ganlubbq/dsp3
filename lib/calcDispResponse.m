% Calculate frequency response of optical chromatic dispersion with
% accumulated dispersion and dispersion slop
% 
% Example: 
%   [H] = calcDispResponse(nSample,fs,lambda,lambda0,DL,SL)
% 
% Input: 
%       nSample     - number of samples
%       fs          - sampling frequency [Hz]
%       lambda      - center wavelength of signal
%       lambda0     - center wavelength of global sim
%       DL          - accumulated linear dispersion in [s/m]
%       SL          - acuumulated dispersion slop in [s2/m]
% 
% Reference: 
% 
% Note: The higher-order effect only needed for very high bandwidth
% 
% See also: calcBesselRespon, calcOptGaussFlt
% 
% Copyright 2015 Dawei Wang 
% 
% Email: dawei.zju@gmail.com 

function [H] = calcDispResponse(nSample,fs,lambda,lambda0,DL,SL)

c = 299792458;

% frequency grid in Hz
freqGrid = [(0:nSample/2-1)'; flipud(-(1:(nSample/2))')] *fs/nSample;

% DL of multiband
DL = DL + (lambda-lambda0)*SL;

dPhase = lambda^2 * DL * pi / c * freqGrid.^2 + ...
    -lambda^3 * DL * 2*pi / 3 / c^2 * freqGrid.^3 + ...
    -lambda^4 * SL * pi / 3 / c^2 * freqGrid.^3;

% frequency response of linear dispersion
H = exp(1j*dPhase);

return


