% Modulate single polarization optical wave using IQ MZM
% 
% Example: oSigMod = oeModIqNested(oSig,eSigI,eSigQ,ERdB,Vpi,V1,V2,V3)
% 
% Input: 
%		V1		- Vpi/2
%		V2		- Vpi/2
%		V3		- Vpi/2
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright 2015 Dawei Wang 
% 
% Email: dawei.zju@gmail.com 

function oSigMod = oeModIqNested(oSig,eSigI,eSigQ,ERdB,Vpi,V1,V2,V3)

oSig = oSig(:);
eSigI = eSigI(:);
eSigQ = eSigQ(:);

Imperfect          = 1/sqrt(10^(ERdB/10));
psr1       = sqrt(0.5+Imperfect);
psr2       = sqrt(0.5-Imperfect);

% push-pull
phi1 = +(pi.*eSigI./Vpi + pi*V1/Vpi);
phi2 = -(pi.*eSigI./Vpi + pi*V1/Vpi);
oSigMod1 = oSig./sqrt(2).*(psr1*exp(1j*phi1) + psr2*exp(1j*phi2) );

phi1 = +(pi.*eSigQ./Vpi + pi*V2/Vpi);
phi2 = -(pi.*eSigQ./Vpi + pi*V2/Vpi);
oSigMod2 = oSig./sqrt(2).*(psr1*exp(1j*phi1) + psr2*exp(1j*phi2) );

oSigMod = oSigMod1 + exp(1j*V3*pi/Vpi)*oSigMod2;

return

