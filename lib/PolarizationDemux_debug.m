function [argout err h1 h2] = PolarizationDemux_debug(signal,mn,sps, ...
    polm, mu, ntaps, errid, iter, appLMS )
%POLARIZATIONDEMUX Polarization demultiplexing and channel equalization
%   algorithm includes constant modulus algorithm (CMA) and least mean
%   squared (LMS). Multi-stage configuration is applied.
%
%   Example
%   
%   See also CMA_FILTER.c

%   Copyright2011 default

%   Last Modified 30/6/2012

x = DspAlg.Normalize(signal, mn);
x = x / (sqrt(mn)-1);
nt = ntaps(1);
% taps initialization
halfnt = floor( nt/2 );
hzero = zeros(nt, 2, 2);
if ~polm
    r = mean( x(:,1) ./ x(:,2) );
    M = rotpolar(1, r).';
else
    M = [1 0; 0 1];
end
hzero( halfnt+1, :, :) = M;
h1 = squeeze( hzero(:, 1, :) );
h2 = squeeze( hzero(:, 2, :) );
if nt == 1
    h1 = h1.';
    h2 = h2.';
end
cstl = get_constellation(mn);
method = 0;
extendx = [ x(end-halfnt+1:end,:); x ; x(1:halfnt,:) ];
% start all 3 stages
for stage = 1:3
    if ntaps(stage)>nt
        temp = (ntaps(stage)-nt)/2;
        h1 = [zeros(temp,2);h1;zeros(temp,2)];
        h2 = [zeros(temp,2);h2;zeros(temp,2)];
        halfnt = floor(ntaps(stage) / 2);
        extendx = [ x(end-halfnt+1:end,:); x ; x(1:halfnt,:) ];
        nt = ntaps(stage);
    end
    cm = get_constantmodulus(cstl, errid(stage));
    while iter(stage)
        [xx h1 h2 err] = DspAlg.CmaEqualizer2_debug_mcma(extendx,h1,h2,nt,mu(stage),cm,sps,errid(stage),stage,method);
        iter(stage) = iter(stage) - 1;
    end
end
% start LMS
if appLMS
    for n = 1:10
        [xx mse deth] = DspAlg.DDLMS_FILTER(extendx,h1,h2,nt,mu(3),cstl,sps);
    end
end
% format the output
argout = xx(1:sps:end,:);
err = err(1:sps:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = get_constellation(mn)
switch mn
    case 1
        c = [complex(-0,0); complex(1,0)];
    case 2
        c = [complex(-1,0); complex(1,0)];
    case 4
        x = [+1 -1]; quar = [x;x]; quar = quar(:);
        y = [-1 +1]; inph = [y y]; inph = inph.';
        c = inph + 1j*quar;
    case 16
        x = [3 1 -1 -3]; quar = [x;x;x;x]; quar = quar(:);
        y = [-3 -1 1 3]; inph = [y y y y]; inph = inph.';
        d = inph + 1j*quar;
        c = d/3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = get_constantmodulus(c, eid)
switch eid
    case 1
        r = mean(abs(c).^4) / mean(abs(c).^2);
    case 2
        r = mean(real(c).^4) / mean(real(c).^2);
    case 3
        r = [sqrt(1+1); sqrt(1+9); sqrt(9+9)]/3;
    case {4,5,6,7}
        r = [sqrt(1+1); sqrt(1+9); sqrt(9+9)].^2/9;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = rotpolar( x , r )
if abs(r) < 0.5
    m = abs(r);
    delta = angle(r);
    alpha = m^2 / (m^2 + 1);
else
    m = abs( 1 / r);
    delta = -angle( 1 / r );
    alpha = 1 / ( m^2 + 1 );
end
M = [sqrt(alpha)*exp(-1j*delta) -sqrt(1-alpha)*exp(-1j*delta);...
     sqrt(1-alpha)               sqrt(alpha) ];
y = x * M;