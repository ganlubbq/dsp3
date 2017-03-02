function h = chessboard(c)

% 1 is white and 0 is black
h = pcolor(c);

% this is default setting
% shading faceted

colormap(gray(2));

axis ij
axis square

return
