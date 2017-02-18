% Calculate Jones vector with parameter of azimuth and ellipticity in
% degree
% 
% Example: jv = calcJonesVector(azi,ell)
% 
% Input: 
% 
% Reference: 
% 
% Note: 
% 
% See Also: 
% 
% Copyright 2015 Default

function jv = calcJonesVector(azi,ell)

while abs(azi) > 90
    azi = (abs(azi)-180) * sign(azi);
end
while abs(ell) > 45
    ell = (abs(ell)-90) * sign(ell);
end

% convert to radius
ita = azi / 180 * pi;
eps = ell / 180 * pi;

% get the power split ratio
k = ( 1 - cos(2*ita)*cos(2*eps) ) / 2;

% get the phase difference
if k == 0 || k == 1
    d = 0;
else
    d = asin( sin(2*eps) / 2 / sqrt(k*(1-k)) );
end

jv = [sqrt(1-k);    sqrt(k) * exp(1j*d)];

% if it is in Western Hemisphere %
if ita<0
    jv = [sqrt(1-k);        -sqrt(k) * exp(1j*-d)];
end

return

