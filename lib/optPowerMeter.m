%OPTPOWERMETER Optical power meter
%   Define,
%   Ps: the power of one sample
%   Ts: the duration of one sample
%   P = E[Ps*Ts]/Ts = E[Ps] = E[abs(a)^2]
% 
% Example: 
% 
% Input: 
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

function ydisp = optPowerMeter( x, dbm)

if nargin<2
    dbm = 0;
end

a = x.E;
y = mean(sum(abs(a).^2));

if dbm
    ydisp = 10 * log10(1000*y);
    % fprintf('%.4f%s\n', ydisp, ' dBm');
else
    ydisp = 1000*y;
    % fprintf('%.4f%s\n', ydisp, ' mw');
end

return

