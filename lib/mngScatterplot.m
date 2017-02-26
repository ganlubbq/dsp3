function mngScatterplot(input1, input2, varargin)
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

h1 = scatterplot(input1, varargin{:}); grid on
h2 = scatterplot(input2, varargin{:}); grid on

% ch = get(get(h1,'Children'),'Children');
% set(ch,'Marker','.','MarkerSize',3) 
% ch = get(get(h2,'Children'),'Children');
% set(ch,'Marker','.','MarkerSize',3) 

% ch = get(h1,'Children');
% set(ch,'Color',[0 0 0])
% ch = get(h2,'Children');
% set(ch,'Color',[0 0 0])

pos1 = get(h1,'Position');
pos1 = pos1 - [pos1(3)/2 0 0 0];
pos2 = pos1 + [pos1(3)+15 0 0 0];
set(h1,'Position',pos1)
set(h2,'Position',pos2)

drawnow;

return
