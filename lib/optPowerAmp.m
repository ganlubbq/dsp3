% DESCRIPTION
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

function [y] = optPowerAmp( x, power )

y = copy(x);

a = optPowerMeter(x);

y.E = x.E.* sqrt(power/a);

return

