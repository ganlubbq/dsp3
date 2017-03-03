function h = chessboard(c)

map_c = [1.0  1.0  1.0;
    0.5  0.5  0.5];

%
h = pcolor(c);

% this is default setting
% shading faceted

% smallest element in c will take the first color and largest take the last
% color
colormap(map_c);

% more intuitive axis orientation
axis ij
% axis square

return
