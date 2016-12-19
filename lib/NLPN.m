function CxIQ = NLPN(CxIQ,Lspan,Att_dB_m,Gamm,Ps,NSpans)
% This procedure demonstrates the algorithm for nonlinear phase noise
%   compensation described in [1]. 
%   1. Keang-Po Ho, Josef M Kahn "Electronic Compensation
%   Technique to Mitigate Nonlinear Phase Noise"

%   length of one fiber span [m] : Lspan
%   fiber attenuation [dB/m] : Att_dB_m
%   nonlinear coefficient 2*pi*NonlinearIndex*ReferenceFrequency/c/CoreArea [1/(W*m)] : Gamm
%   received power W : Ps
%   number of fiber spans : NSpans

%   fiber Loss
Att = Att_dB_m/4.343;
%   effective fiber length
if Att>0
    Leff = (1-exp(-Att.*Lspan))./Att;
else
    Leff = Lspan;
end;

%   correction factor
alpha = Gamm*Leff*(NSpans+1)/2;
%   power scaling factor
kP = Ps/((CxIQ*CxIQ')/length(CxIQ));

%   perform correction
CxIQ = CxIQ.*exp(-i*CxIQ.*conj(CxIQ).*kP.*alpha);