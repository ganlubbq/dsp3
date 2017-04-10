function [y V] = endlessPolRot(x)
%ENDLESSPOLROT Endless polarization rotation of a optical signalDawei Wang 
% 
%   The rotation speed is denoted by V with unit of [rad/s] and Fs is the
%   sampling frequency
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
% V = V - mod(2*pi/x.T,V);
V = 2*pi/x.T/1;

y = copy(x);

delta = V/x.fs;
L = length(x.E);
theta = (1:L)*delta;

Ox = x.Ex.*cos(theta) - x.Ey.*sin(theta);
Oy = x.Ex.*sin(theta) + x.Ey.*cos(theta);

y.E = [Ox;Oy];

return
