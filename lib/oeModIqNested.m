function oSigMod = oeModIqNested(oSig, eSigI, eSigQ, ERdB, Vpi, V1, V2, V3)
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


oSig = oSig(:);
eSigI = eSigI(:);
eSigQ = eSigQ(:);

Imperfect  = 1 / sqrt(10 .^ (ERdB / 10));
psr1       = sqrt(0.5 + Imperfect);
psr2       = sqrt(0.5 - Imperfect);

% push-pull
phi1 =   pi * (eSigI ./ Vpi + V1 ./ Vpi);
phi2 = - pi * (eSigI ./ Vpi + V1 ./ Vpi);
oSigMod1 = oSig ./ sqrt(2) .* (psr1 * exp(1i * phi1) + psr2 * exp(1i * phi2));

phi1 =   pi * (eSigQ ./ Vpi + V2 ./ Vpi);
phi2 = - pi * (eSigQ ./ Vpi + V2 ./ Vpi);
oSigMod2 = oSig ./ sqrt(2) .* (psr1 * exp(1i * phi1) + psr2 * exp(1i * phi2));

oSigMod = oSigMod1 + exp(1i * V3 * pi / Vpi) * oSigMod2;

return
