function mngFigureWindow( h1, h2, config )
% Manage position of 2 figures in the mode configured by input
% 
% Example: mngFigureWindow( h1, h2, config )
% 
% Input: config - lr or ud
% 
% Reference: 
% 
% Note: 
% 
% See Also: 

if nargin < 3
    config = 'lr';
end

scrsz = get(0,'ScreenSize');

p1 = get(h1,'Position');
p2 = get(h2,'Position');

if strcmpi(config,'lr')
    p1(1) = 100;
    p2(1) = p1(1) + p1(3) + 30;
    set(h1,'Position',p1)
    set(h2,'Position',p2)
elseif strcmpi(config,'ud')
    p1(2) = scrsz(4) - p1(4) - 80;
    p2(2) = p1(2) - p2(4) - 80;
    set(h1,'Position',p1)
    set(h2,'Position',p2)
end

drawnow;

return

