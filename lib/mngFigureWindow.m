function mngFigureWindow(h1, h2)
% Manage position of 2 figures in the mode configured left-right side by
% side
% 
% Example: mngFigureWindow(h1, h2, config)
% 
% See Also: 
scrsz = get(0,'ScreenSize');

screen_width = scrsz(3);
% screen_height = scrsz(4);

p1 = get(h1,'Position');
p2 = get(h2,'Position');

p1(1) = round(screen_width/2) - p1(3);
p2(1) = p1(1) + p1(3) + 10;
set(h1,'Position',p1);
set(h2,'Position',p2);

drawnow;

return
