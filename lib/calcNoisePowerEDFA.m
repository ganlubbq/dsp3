function [psd] = calcNoisePowerEDFA(NA,gain,lambda,nsp)
% Calculate noise power density per state of polarization for periodic EDFA
% link
% 
% Example: [psd] = calcNoisePowerEDFA(NA,gain,lambda,nsp)
% 
% Input: 
%       NA          - number of EDFA stage
%       gain        - linear gain
%       lambda      - center wavelength of signal
%       nsp         - noise spontaneous emission factor >=1
% 
% Reference: Essiambre, et al "Capacity Limits of Optical Fiber Networks,"
% in Lightwave Technology, Journal of , vol.28, no.4, pp.662-701, Feb.15,
% 2010
% 
% Note: EDFA noise figure F can be expressed as F=2*nsp-(2*nsp-1)/G which
% is about 2*nsp if the gain G is very large
% 
% See Also: calcNoisePowerIDRA
% 
% Copyright 2015 Default

% Planck constant [J*s]
h = 6.626068e-34;

c = 299792458;

v = c/lambda;

% noise power spectrum density [W/Hz]
psd = NA*(gain-1)*h*v*nsp;

return
