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
% Copyright 2015 default

function [y] = optPowerAmp( x, power )

y = copy(x);

a = optPowerMeter(x);

y.E = x.E.* sqrt(power/a);

return

