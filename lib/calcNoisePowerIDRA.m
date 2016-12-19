% Calculate noise power density per state of polarization for ideal Raman
% amplifier
% 
% Example: [psd] = calcNoisePowerIDRA(alpha,totlen,lambda,KT)
% 
% Input: 
%       alpha       - linear loss factor
%       totlen      - total fiber length [m]
%       lambda      - center wavelength of signal
%       KT          - photon occupancy factor = 1.13 at 1550nm and room
% 
% Reference: Essiambre, et al "Capacity Limits of Optical Fiber Networks,"
% in Lightwave Technology, Journal of , vol.28, no.4, pp.662-701, Feb.15,
% 2010
% 
% Note: 
% 
% See Also: calcNoisePowerEDFA
% 
% Copyright 2015 Dawei Wang 
% 
% Email: dawei.zju@gmail.com 

function [psd] = calcNoisePowerIDRA(alpha,totlen,lambda,KT)

% Planck constant [J*s]
h = 6.626068e-34;

c = 299792458;

v = c/lambda;

% noise power spectrum density [W/Hz]
psd = alpha*totlen*h*v*KT;

return

